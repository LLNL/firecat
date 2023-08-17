from firecat import BicarbGMPNPSolver, potential_continuation
from echemfem import IntervalBoundaryLayerMesh
import argparse
import os
import numpy as np
from petsc4py import PETSc

"""
This code pproduces the 1D results in:
N. Govindarajan, T. Y. Lin, T. Roy, C. Hahn, J. B. Varley. Coupling
microkinetics with continuum transport models to understand electrochemical CO2
reduction in flow reactors. To appear in PRX Energy, 2023.
"""

pprint = PETSc.Sys.Print
parser = argparse.ArgumentParser()
parser.add_argument('-results', type=str, default='results')
parser.add_argument('-U_final', type=float, default=1.15)
parser.add_argument('-tol', type=float, default=1e-3)
parser.add_argument('-bdlayer', type=float, default=64.7)
args, unknown = parser.parse_known_args()
pprint(args)


# set up boundary conditions
electrode_id = 2
boundary_markers = {"bulk dirichlet": (1,),  #C = C_0
		     "bulk": (1,), # U_liquid = 0
		     "neumann": (electrode_id,), 
		     "robin": (electrode_id,),
		     }

delta = args.bdlayer * 1e-6
mesh = IntervalBoundaryLayerMesh(200, delta, 800, 1e-8)

# set up continuum model
solver = BicarbGMPNPSolver(mesh, boundary_markers, p=2)
solver.save_solutions=False
solver.U_app.assign(0)
solver.setup_solver()
solver.electrode_id = electrode_id
solver.solve()

# Defining the mkm model
args.mkm_model = 'catmap_CO2R'

# array of potentials to consider
Us = np.linspace(-0.3, -0.6, 4)
Us = np.concatenate((Us, np.linspace(-0.65, -1.0, 8)))
Us = np.concatenate((Us, np.linspace(-1.01, -1.15, 15)))
#args.theta = [0.5, 0.25]

# coupling with catmap
potential_continuation(solver, Us, args)
