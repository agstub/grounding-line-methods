# This file contains the functions needed for solving the Stokes system.
from params import rho_i,g,tol,B,rm2,rho_w,C,eps_p,eps_v,sea_level,dt,quad_degree,Lngth
from boundaryconds import mark_boundary,apply_bcs
from hydrology import Vdot, sl_change
import numpy as np
from dolfin import *

def dPi(u,nu):
        # Derivative of penalty functional for enforcing impenetrability on the ice-bed boundary.
        un = dot(u,nu)
        return un+abs(un)

def Pi(u,nu):
        # Penalty functional for enforcing impenetrability on the ice-bed boundary.
        un = dot(u,nu)
        return 0.5*(un**2.0+un*abs(un))


def weak_form_lake(u,p,pw,v,q,qw,f,g_lake,ds,nu,T,lake_vol_0,t):
    # Weak form of the subglacial lake problem

    # Measures of the lower boundary and ice-water boundary
    L0 = Constant(assemble(1*ds(4))+assemble(1*ds(3))) # entire lower boundary
    L1 = Constant(assemble(1*ds(4)))                   # ice-water boundary

    # Nonlinear residual
    Fw = B*((inner(sym(grad(u)),sym(grad(u)))+Constant(eps_v))**(rm2/2.0))*inner(sym(grad(u)),sym(grad(v)))*dx\
         +(- div(v)*p + q*div(u))*dx - inner(f, v)*dx\
         + (g_lake+pw+Constant(rho_w*g*dt)*(dot(u,nu)+Constant(Vdot(lake_vol_0,t)/L1)))*inner(nu, v)*ds(4)\
         + qw*(inner(u,nu)+Constant(Vdot(lake_vol_0,t))/(L0))*ds(4)\
         + (g_lake+pw+Constant(rho_w*g*dt)*(dot(u,nu)+Constant(Vdot(lake_vol_0,t)/L1)))*inner(nu, v)*ds(3)\
         + qw*(inner(u,nu)+Constant(Vdot(lake_vol_0,t))/(L0) )*ds(3)\
         + Constant(1/eps_p)*dPi(u,nu)*dot(v,nu)*ds(3)\
         + Constant(C)*((Constant(eps_v)+inner(dot(T,u),dot(T,u)))**(rm2/2.0))*inner(dot(T,u),dot(T,v))*ds(3)
    return Fw


def stokes_solve_lake(mesh,lake_vol_0,s_mean,t):
        # Stokes solver using Taylor-Hood elements and a Lagrange multiplier
        # for the water pressure.

        # Define function spaces
        P1 = FiniteElement('P',triangle,1)     # Pressure
        P2 = FiniteElement('P',triangle,2)     # Velocity
        R  = FiniteElement("R", triangle,0)    # Mean water pressure
        element = MixedElement([[P2,P2],P1,R])
        W = FunctionSpace(mesh,element)        # Function space for (u,p,pw)

        #---------------------Define variational problem------------------------
        w = Function(W)
        (u,p,pw) = split(w)             # (velocity,pressure,mean water pressure)
        (v,q,qw) = TestFunctions(W)     # test functions corresponding to (u,p,pw)

        # Define Neumann condition at ice-water interface
        g_lake = Expression('rho_w*g*(s_mean-x[1])',rho_w=rho_w,g=g,s_mean=s_mean,degree=1)

        f = Constant((0,-rho_i*g))        # Body force
        nu = FacetNormal(mesh)            # Outward-pointing unit normal to the boundary
        I = Identity(2)                   # Identity tensor
        T = I - outer(nu,nu)              # Orthogonal projection (onto boundary)

        # Mark the boundary and define a measure for integration
        boundary_markers = mark_boundary(mesh)
        ds = Measure('ds', domain=mesh, subdomain_data=boundary_markers)

        # Define weak form
        Fw = weak_form_lake(u,p,pw,v,q,qw,f,g_lake,ds,nu,T,lake_vol_0,t)

        #Apply Dirichlet BC on side walls
        bcs_u =  apply_bcs(W,boundary_markers)

        # Solve for (u,p,pw).
        solve(Fw == 0, w, bcs=bcs_u,solver_parameters={"newton_solver":{"relative_tolerance": 1e-14,"maximum_iterations":100}},form_compiler_parameters={"quadrature_degree":quad_degree,"optimize":True,"eliminate_zeros":False})

        # Compute penalty functional residiual
        P_res = assemble(Pi(u,nu)*ds(3))

        # Return solution w and penalty functional residual Perr
        return w,P_res

def weak_form_marine(u,p,v,q,f,g_base,g_out,ds,nu,T):
    # Weak form of the marine ice sheet problem

    Fw = B*((inner(sym(grad(u)),sym(grad(u)))+Constant(eps_v))**(rm2/2.0))*inner(sym(grad(u)),sym(grad(v)))*dx\
         + (- div(v)*p + q*div(u))*dx - inner(f, v)*dx\
         + (g_base+Constant(rho_w*g*dt)*inner(u,nu))*inner(nu, v)*ds(4)\
         + Constant(C)*((inner(dot(T,u),dot(T,u))+Constant(eps_v))**(rm2/2.0))*inner(dot(T,u),dot(T,v))*ds(3)\
         + Constant(1.0/eps_p)*dPi(u,nu)*dot(v,nu)*ds(3)\
         + (g_base+Constant(rho_w*g*dt)*inner(u,nu))*inner(nu, v)*ds(3)\
         + g_out*inner(nu, v)*ds(2)
    return Fw


def stokes_solve_marine(mesh,F_h,t):
        # Stokes solver using Taylor-Hood elements.

        # Define function spaces
        P1 = FiniteElement('P',triangle,1)     # Pressure
        P2 = FiniteElement('P',triangle,2)     # Velocity
        element = MixedElement([[P2,P2],P1])
        W = FunctionSpace(mesh,element)        # Function space for (u,p)

        #---------------------Define variational problem------------------------
        w = Function(W)
        (u,p) = split(w)
        (v,q) = TestFunctions(W)

        # Neumann condition at outflow boundary
        h_out = float(F_h(Lngth))        # Surface elevation at outflow boundary
        g_out = Expression('rho_i*g*(h_out-x[1])',rho_i=rho_i,g=g,h_out=h_out,degree=1)

        # Neumann condition at ice-water boundary
        g_base = Expression('rho_w*g*(sea_level-x[1])',rho_w=rho_w,g=g,sea_level=sea_level+sl_change(t),degree=1)

        f = Constant((0,-rho_i*g))        # Body force
        nu = FacetNormal(mesh)            # Outward-pointing unit normal to the boundary
        I = Identity(2)                   # Identity tensor
        T = I - outer(nu,nu)              # Tangential projection operator

        # Mark bounadries of mesh and define a measure for integration
        boundary_markers= mark_boundary(mesh)
        ds = Measure('ds', domain=mesh, subdomain_data=boundary_markers)

        # Define weak form and apply boundary conditions on the inflow boundary

        bcs_u =  apply_bcs(W,boundary_markers)    # Apply Dirichlet BC


        # Solve for (u,p).
        Fw = weak_form_marine(u,p,v,q,f,g_base,g_out,ds,nu,T)


        solve(Fw == 0, w, bcs=bcs_u,solver_parameters={"newton_solver":{"relative_tolerance": 1e-14,"maximum_iterations":100}},form_compiler_parameters={"quadrature_degree":quad_degree,"optimize":True,"eliminate_zeros":False})

        # Compute penalty functional residiual
        P_res = assemble(Pi(u,nu)*ds(3))

        return w,P_res


def get_zero_m(mesh):
        # Get zero element of function space for marine ice sheet problem.
        # Only used for setting initial conditions; see main.py.

        # Define function spaces
        P1 = FiniteElement('P',triangle,1)     # Pressure
        P2 = FiniteElement('P',triangle,2)     # Velocity
        element = MixedElement([[P2,P2],P1])
        W = FunctionSpace(mesh,element)        # Function space for (u,p)

        w = Function(W)

        return w


def get_zero_l(mesh):
        # Get zero element of function space for subglacial lake problem.
        # Only used for setting initial conditions; see main.py.

        # Define function spaces
        P1 = FiniteElement('P',triangle,1)     # Pressure
        P2 = FiniteElement('P',triangle,2)     # Velocity
        R = FiniteElement("R", triangle,0)     # Mean water pressure
        element = MixedElement([[P2,P2],P1,R])
        W = FunctionSpace(mesh,element)        # Function space for (u,p,pw)
        w = Function(W)

        return w
