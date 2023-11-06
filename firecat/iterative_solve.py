from firedrake import *
from firecat import run_catmap, convert_TOF_to_flux, __file__
import os
import numpy as np
from petsc4py import PETSc
pprint = PETSc.Sys.Print
## Everything is hard-coded for gold for now

def get_activities(solver, args):
    """ Returns a list of activities for CO2 and H 
    solver: EchemFEM solver
    args: object containing various properties such as
            dofs:       indices of the dofs of electrode
    """
    dofs = args.dofs
    cCO2, cOH, cHCO3, cCO32, cH, cK, Phi = solver.u.subfunctions
    aCO2 = cCO2.vector().get_local()[dofs[0]] / solver.conc_params[0]["bulk"]
    aH = cH.vector().get_local()[dofs[0]] / 1000.
    return [aCO2, aH]

def prate_to_fluxes(prate, local_flux, local_fluxOH, dof):
    """ takes prate and updates fluxes
    prate: production rate
    local_flux: Local data for CO2 flux (numpy.array)
    local_fluxOH: Local data for OH flux (numpy.array)
    dof: index of the dof
    """
    tof = -prate
    flux = convert_TOF_to_flux(tof) * 1e4 # mol/m2
    local_flux[dof] = flux
    local_fluxOH[dof] = -2 * flux

def iterative_coupling(solver, args, U, Upzc, u_old0, u_old1):
    """ For a given potential, solve the echemfem - catmap coupled problem
    solver: EchemFEM solver
    U: potential value
    args: object containing various properties such as
            mkm_model:  name of mkm model
            dofs:       indices of the dofs of electrode
            tol:        relative tolerance for the iterative coupling
    Upzc: potential of zero charge
    u_old0: EchemFEM solution of the previous iteration
    u_old1: EchemFEM solution from two iterations before
    """
    if getattr(args, 'template_root', None):
        template_root = args.template_root
    else:
        template_root = os.path.dirname(__file__) + '/mkm_models/'

    flux_diff_rel = 1e10
    i = 0
    if getattr(args, 'theta', None):
        theta0_ = args.theta[0]
        theta1_ = args.theta[1]
    else:
        theta0_ = 0.35
        theta1_ = 0.1
    dofs = args.dofs

    C_gap = solver.physical_params.get("gap capacitance constant")
    if C_gap is None:
        C_gap = solver.physical_params["gap capacitance"]

    #local copy of wall fluxes
    flux = Function(solver.V)
    fluxOH = Function(solver.V)
    flux_old = Function(solver.V)


    while flux_diff_rel > args.tol:
        pprint("Voltage = " + str(np.round(U,8)) +" V")
        pprint(i)
        i += 1
        pprint("Number of iterations: ", i)
        pprint(flux_diff_rel)
        if i > 10 and i < 40:
            theta0 = 0.55
            theta1 = 0.3
        elif i >= 40:
            theta0 = 0.6
            theta1 = 0.35
        else:
            theta0 = theta0_
            theta1 = theta1_ 

        solver.u.assign(solver.u * (1 - theta0 - theta1) + u_old0 * theta0 + u_old1 * theta1)
        activities = get_activities(solver, args) # contains aCO/aCO2 and aOH

        while any((a<0).any() for a in activities):
            pprint("Negative concentration. Increasing damping")
            theta0 += 0.05
            theta1 += 0.05
            if theta0 + theta1 > 0.999:
                pprint("WARNING: OVERDAMPING")
                exit()    
            solver.u.assign(solver.u * (1 - theta0 - theta1) + u_old0 * theta0 + u_old1 * theta1)
            activities = get_activities(solver, args)

        aH = activities[1] # hard-coded for now
        pH = -np.log10(aH)

        if args.mkm_model == 'catmap_CO2R':
            _, _, _, _, _, _, Phi = solver.u.subfunctions
            sigma = 1e2*C_gap * (U - Upzc - Phi.vector().get_local()[dofs[0]]) #muC/cm2
            surfacePhi = Phi.vector().get_local()[dofs[0]]
        else:
            sigma=None
            surfacePhi=None

        flux_old.assign(flux)
        flux.assign(Constant(0.))
        fluxOH.assign(Constant(0.))

        local_flux = flux.vector().get_local() 
        local_fluxOH = flux.vector().get_local()
        for j in range(dofs[0].size):
            prate = run_catmap(U, pH[j], activities[0][j], sigma=sigma[j], phi=surfacePhi[j], j=j, template_root=template_root)
            prate_to_fluxes(prate, local_flux, local_fluxOH, dofs[0][j])
        local_flux = flux.vector().set_local(local_flux) 
        local_fluxOH = flux.vector().set_local(local_fluxOH)

        if args.mkm_model == 'catmap_CO2R':   
            solver.fluxCO2.assign(flux)
        else:
            solver.fluxCO.assign(flux)
        solver.fluxOH.assign(fluxOH)
        
        if hasattr(solver, "bump"):
            flux_diff_rel = norm((flux_old - flux)*solver.bump)/norm(flux*solver.bump)
        else:
            flux_diff_rel = norm(flux_old - flux)/norm(flux)

        u_old1.assign(u_old0)
        u_old0.assign(solver.u)
        solver.solve()

    solver.output_state(solver.u, prefix=args.results+'/'+str(np.round(U,8))+'/')
    if getattr(args, 'save_checkpoints', None):
        with CheckpointFile("checkpoint_" + str(np.round(U,8)) + ".h5", 'w') as afile:
            afile.save_mesh(solver.mesh)
            afile.save_function(u_old0)
            afile.save_function(solver.u)

    return u_old0, u_old1

def potential_continuation(solver, Us, args):
    """ Solve the coupled problem for a list of potentials
    solver: EchemFEM solver
    Us: array of applied potential values
    args: object containing various properties such as
            mkm_model:  name of mkm model
            dofs:       indices of the dofs of electrode
            tol:        relative tolerance for the iterative coupling
    """

    if solver.physical_params.get("Upzc"):
        Upzc = solver.physical_params.get("Upzc")
    else: # if using Robin BC
        Upzc = 0.

    # Only works for CG function spaces
    def get_boundary_dofs(V, i):
        u = Function(V)
        bc = DirichletBC(V, Constant(1), i)
        bc.apply(u)
        return np.where(u.vector().get_local()[:]==1)

    # get the indices of the boundary points
    args.dofs = get_boundary_dofs(solver.V, solver.electrode_id)
    if hasattr(solver, "bump"):
        np.save('testing/bump_values.npy', solver.bump.vector()[args.dofs])

    # array of potentials to consider
    if hasattr(args, "U_final"):
        Us = Us[Us>=-args.U_final]
        if -args.U_final not in Us: Us = np.append(Us,-args.U_final)
        if getattr(args, 'checkpoint', None):
            Us = Us[Us<-args.checkpoint-1e-16]

    pprint("All potentials = ", Us)

    fluxes = []

    j_COs = []
    F = 96485

    # initialize the previous iterations
    if getattr(args, 'checkpoint', None) is None:
        solver.U_app.assign(Us[0] - Upzc)
        solver.solve()
        u_old0 = Function(solver.W, name="u_old0").assign(solver.u)
        u_old1 = Function(solver.W, name="u_old1").assign(solver.u)
    else:
        u_old0 = args.u_old0
        u_old1 = args.u_old1
        solver.u.assign(u_old0)
    # loop through potentials
    for U in Us:

        solver.U_app.assign(U - Upzc)

        u_old0, u_old1 = iterative_coupling(solver, args, U, Upzc, u_old0, u_old1)
            
        cCO2, cOH, cHCO3, cCO32, cH, cK, Phi = solver.u.subfunctions

        j_CO = -2 * solver.neumann(cCO2, {"name": "CO2"}, solver.u) * F * 1e-4 * 1e3   # mA/cm2
        j_CO_avg = assemble(j_CO * ds(solver.electrode_id)) \
                    / assemble(Constant(1) * ds(solver.electrode_id, domain=solver.mesh))
        j_COs.append(j_CO_avg)

        pprint("j_CO_avg: ",j_CO_avg)

        if solver.mesh.comm.rank == 0:
            if not os.path.exists('testing'):
               os.makedirs('testing')
            filename = "testing/jCO.txt"
            if U == Us[0] and getattr(args, 'checkpoint', None) is None:
                file_ = open(filename, "w")
            else:
                file_ = open(filename, "a")
            print("U = ", U,", j_CO = ", j_CO_avg, file=file_)
            file_.close()

        File("testing/jCO.pvd").write(Function(solver.V, name='jCO').interpolate(j_CO))
