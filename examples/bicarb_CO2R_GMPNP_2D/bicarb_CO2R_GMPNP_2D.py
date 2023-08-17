from firecat import BicarbGMPNPSolver, potential_continuation
from firedrake import Mesh, SpatialCoordinate, as_vector, Constant
import argparse
import os
import numpy as np
from petsc4py import PETSc

"""
This code produces the 2D results in:
N. Govindarajan, T. Y. Lin, T. Roy, C. Hahn, J. B. Varley. Coupling
microkinetics with continuum transport models to understand electrochemical CO2
reduction in flow reactors. To appear in PRX Energy, 2023.
"""

pprint = PETSc.Sys.Print
parser = argparse.ArgumentParser()
parser.add_argument('-mesh', type=str, default='firecat2D_1cm_1.msh')
parser.add_argument('-results', type=str, default='results')
parser.add_argument('-U_final', type=float, default=1.1)
parser.add_argument('-tol', type=float, default=1e-3)
parser.add_argument('-shear_rate', type=float, default=1.91)
parser.add_argument('-height_factor', type=float, default=None)
parser.add_argument('-checkpoint', type=float, default=None)
parser.add_argument('-save_checkpoints', type=int, default=1)
args, unknown = parser.parse_known_args()
pprint(args)

if args.checkpoint is not None:
    chk_name = "checkpoint_-" + str(args.checkpoint) + ".h5"
    pprint("loading file " + chk_name)
    with CheckpointFile(chk_name, 'r') as afile:
        mesh_chk = afile.load_mesh(os.path.basename(args.mesh))
        u_old0 = afile.load_function(mesh_chk, "solution")
        u_old1 = afile.load_function(mesh_chk, "u_old0") 
        u_old0.rename("u_old0")
        u_old1.rename("u_old1")
        args.u_old0 = u_old0
        args.u_old1 = u_old1

if args.checkpoint is None:
    mesh = Mesh(args.mesh, name=args.mesh)
    if args.height_factor is not None:
        mesh.coordinates.dat.data[:, 1] *= args.height_factor          
else:
    mesh = mesh_chk

# set up boundary conditions
electrode_id = 4
boundary_markers = {"inlet": (1,), #x=0
                         "bulk": (3,),
                         "outlet": (2,), #x=Lx
                         "neumann": (electrode_id,), #y=0
                         "robin": (electrode_id,), #y=0
                         }                           

# set up shear flow velocity field
_, y = SpatialCoordinate(mesh)
velocity = as_vector([args.shear_rate*y, Constant(0)]) # m/s

# set up continuum model
solver = BicarbGMPNPSolver(mesh, boundary_markers, velocity=velocity)
solver.save_solutions=False
solver.setup_solver()
solver.electrode_id = electrode_id

# Defining the mkm model
args.mkm_model = 'catmap_CO2R'

# array of potentials to consider
Us = np.linspace(-0.4, -0.80, 9)
#Us = np.linspace(-0.4, -0.80, 21)
Us = np.concatenate((Us, np.linspace(-0.81, -0.84, 4)))
Us = np.concatenate((Us, np.linspace(-0.845, -1.0, 32)))
Us = np.concatenate((Us, np.linspace(-1.005, -1.1, 20)))
Us = np.concatenate((Us, np.linspace(-1.105, -1.2, 20)))
Us = np.concatenate((Us, np.linspace(-1.205, -1.3, 20)))

# coupling with catmap
potential_continuation(solver, Us, args)
