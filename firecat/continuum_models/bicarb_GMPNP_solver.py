from firedrake import *
from echemfem import EchemSolver
from petsc4py import PETSc
pprint = PETSc.Sys.Print

"""
Continuum-scale model with bicarbonate bulk reactions
in a solution of KHCO3, with GMPNP for the transport.
"""

KW = (1E-14)*pow(1000,2)   #(mol/m^3)^2
K1 = pow(10,-6.37)*1000      #(mol/m^3)
K2 = pow(10,-10.32)*1000   #(mol/m^3)
K3 = K1/KW              #(mol/m^3)^-2
K4 = K2/KW              #(mol/m^3)^-1
k1f = 8.42e3 / 1000     #m3/mol.s
k5f = 2.3e10 / 1000     #m3/mol.s (Note - typo in Schulz)
k1r = k1f/K3            #mol/(m^3*s)
k2f = 3.71e-2           #1/s
k2r = k2f/K1            #m3/mol.s
k3r = 6.0e9 / 1000      #m3/mol.s
k3f = k3r/K4            #1/s
k4r = 59.44             #1/s
k4f = k4r/K2            #m3/mol.s
k5r = KW*k5f            #mol/m3.s


C_1_inf = 0.034*1000
C_K = .5*1000
Keq = k1f*k3f/(k1r*k3r)

C_3_inf = -C_1_inf*Keq/4 + C_1_inf*Keq/4*sqrt(1+8*C_K/(C_1_inf*Keq))
C_4_inf = (C_K - C_3_inf)/2
C_2_inf = C_3_inf/C_1_inf*k1r/k1f
C_5_inf = (k5r/k5f)/C_2_inf

C_CO2_bulk = C_1_inf
C_OH_bulk = C_2_inf
C_HCO3_bulk = C_3_inf
C_CO32_bulk = C_4_inf
C_H_bulk = C_5_inf
C_K_bulk = C_K

#Make sure that electroneutrality holds by manually adjusting K concentration
netcharge = C_K_bulk+ C_H_bulk -2.*C_CO32_bulk - C_HCO3_bulk - C_OH_bulk
pprint("Total charge is:",netcharge )
C_K_bulk = C_K_bulk - netcharge
netcharge = C_K_bulk+ C_H_bulk -2.*C_CO32_bulk - C_HCO3_bulk - C_OH_bulk
pprint("Total charge is:",netcharge )

def bulk_reaction(y):
    yCO2=y[0]
    yOH=y[1]
    yHCO3=y[2]
    yCO3=y[3]
    yH=y[4]
    
    
    dCO2 = -(k1f)   *yCO2*yOH \
            +(k1r)    *yHCO3 \
            -(k2f)  *yCO2 \
            +(k2r)   *yHCO3*yH

    dOH = -(k1f)       *yCO2*yOH \
                   +(k1r)        *yHCO3 \
                   +(k3f)        *yCO3 \
                   -(k3r)       *yOH*yHCO3\
                   -(k5f)      *yOH*yH\
                   +(k5r)

    dHCO3 = (k1f)        *yCO2*yOH\
                   -(k1r)       *yHCO3\
                   +(k3f)        *yCO3\
                   -(k3r)       *yOH*yHCO3\
                   +(k2f)       *yCO2 \
                   -(k2r)      *yHCO3*yH \
                   +(k4f)       *yCO3*yH\
                   -(k4r)      *yHCO3

    dCO3 = -(k3f) *yCO3 \
                   +(k3r)  *yOH*yHCO3\
                   -(k4f)*yCO3*yH\
                   +(k4r) *yHCO3

    dH = (k2f)   *yCO2 \
                   -(k2r)  *yHCO3*yH \
                   -(k4f)  *yCO3*yH\
                   +(k4r)   *yHCO3\
                   -(k5f)  *yOH*yH\
                   +(k5r)
       
    return [dCO2, dOH, dHCO3, dCO3, dH, 0.]

class BicarbGMPNPSolver(EchemSolver):
    def __init__(self, mesh, boundary_markers, velocity=None, SUPG=False, p=1, family="CG", cylindrical=False):
        self.boundary_markers = boundary_markers
        self.vel = velocity
        
        conc_params = []

        conc_params.append({"name": "CO2",
                            "diffusion coefficient": 1.91E-9,  # m^2/s
                            "bulk": C_CO2_bulk,  # mol/m3
                            "z": 0,
                            "solvated diameter": 0.0
                            })

        conc_params.append({"name": "OH",
                            "diffusion coefficient": 5.29E-9,  # m^2/s
                            "bulk": C_OH_bulk,  # mol/m3
                            "z": -1,
                            "solvated diameter": 0.0
                            })

        conc_params.append({"name": "HCO3",
                            "diffusion coefficient": 1.185E-9,  # m^2/s
                            "bulk": C_HCO3_bulk,  # mol/m3
                            "z": -1,
			    "solvated diameter": 0.0
                            })

        conc_params.append({"name": "CO3",
                            "diffusion coefficient": .92E-9,  # m^2/s
                            "bulk": C_CO32_bulk,  # mol/m3
                            "z": -2,
                            "solvated diameter": 0.0
                            })

        conc_params.append({"name": "H",
                            "diffusion coefficient": 9.311E-9,  # m^2/s
                            "bulk": C_H_bulk,  # mol/m3
                            "z": 1,
                            "solvated diameter": 0.0
                            })

        conc_params.append({"name": "K",
                           "diffusion coefficient": 1.96E-9,  # m^2/s
                           "bulk": C_K_bulk,  # mol/m3
                           "z": 1,
                           "solvated diameter": 8.2e-10 # m
                           })

        flow = ["diffusion", "migration", "poisson", "finite size"]
        if velocity:
            flow.append("advection")
            SUPG=True

        physical_params = {"flow": flow,
                           "F": 96485.,  # C/mol
                           "R": 8.3144598,  # J/K/mol
                           "T": 273.15 + 25.,  # K
                           "bulk reaction": bulk_reaction,
                           "vacuum permittivity": 8.8541878128e-12,  # F/m
                           "relative permittivity": 78.4,
                           "Avogadro constant": 6.02214076e23, #1/mol
                           "U_app": Constant(0), #V
                           "gap capacitance": 0.2, # F/m^2
                           "Upzc": 0.16, #V
                           }

        V = FunctionSpace(mesh, 'CG', p)
        self.fluxCO2 = Function(V)
        self.fluxOH = Function(V)

        super().__init__(conc_params, physical_params, mesh, family='CG', p=p, SUPG=SUPG, cylindrical=cylindrical)
        
        snes_newtonls = {"snes_type": "newtonls",
                "snes_linesearch_type": "l2", 
                "snes_monitor": None,
                "snes_converged_reason": None,
                "snes_rtol": 1e-16,
                "snes_atol": 1e-16,
                "snes_stol": 1e-12, #new addition
                "snes_divergence_tolerance": -1, #uncommented
                "snes_max_it": 50,
                "mat_type": "aij",
                "ksp_type": "preonly",
                "pc_type": "lu",
                "pc_factor_mat_solver_type": "mumps",
                }
        self.init_solver_parameters(custom_solver=snes_newtonls)

    def neumann(self, C, conc_params, u):
        name = conc_params["name"]
        if name in ["HCO3", "CO3", "H", "K"]:
            return Constant(0)
        if name == "CO2":
            return self.fluxCO2
        if name == "OH":
            return self.fluxOH

    # These properties are already defined:
    def set_boundary_markers(self):
        pass

    def set_velocity(self):
        pass
