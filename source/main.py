#-------------------------------------------------------------------------------
# This program solves a nonlinear Stokes flow problem describing groundning line
# dynamics in two settings:
#
# (1) A marine ice sheet, OR ...
# (2) A subglacial lake that is (potentially) undergoing water volume changes.
#
# Choose model setup (1) or (2) in the params.py file.
#
# See the JFM manuscript for a complete formulation of the problem.
#-------------------------------------------------------------------------------
import sys
sys.path.insert(0, './scripts')

from dolfin import *
import matplotlib.pyplot as plt
import numpy as np
from stokes import stokes_solve_lake,stokes_solve_marine,get_zero_m,get_zero_l
from geometry import interface,bed
from meshfcns import mesh_routine
from plotting import *
import scipy.integrate as scpint
import os
from params import (rho_i,g,tol,t_final,nt_per_year,Lngth,Hght,nt,dt,model,
                    print_convergence,X_fine,nx,tides,DX_s,model_setup)

#--------------------Initial conditions-----------------------------------------
# Compute initial mean elevation of ice-water interface and initial lake volume.
s_mean0 = np.mean(interface(X_fine)[interface(X_fine)-bed(X_fine)>tol])
lake_vol_0 = scpint.quad(lambda x: interface(x)-bed(x),0,Lngth,full_output=1)[0]
#-------------------------------------------------------------------------------

os.mkdir('results')   # Make a directory for the results.

if print_convergence == 'off':
    set_log_level(40)    # Suppress Newton convergence information if desired.

# Create VTK files
vtkfile_u = File('results/stokes/u.pvd')
vtkfile_p = File('results/stokes/p.pvd')

# Load mesh
if model == 'marine' and tides=='on':
    meshname = 'tides'+'_DX'+str(int(DX_s))+'.xml'
else:
    meshname = str(model)+'_DX'+str(int(DX_s))+'.xml'

if model_setup == 'wedge_test':
    meshname = 'wedge.xml'

if model == 'marine' and tides=='off':
    # Create initial mesh for tides simulation by running the marine model with
    # tides turned off.
    new_mesh = File('tides_DX+'+str(int(DX_s))+'.xml')

mesh = Mesh('./meshes/'+meshname)


# Define arrays for saving surfaces, lake volume, water pressure, and
# grounding line positions over time.
Gamma_s = np.zeros((nx,nt))       # Basal surface
Gamma_h = np.zeros((nx,nt))       # Upper surface
s_mean = np.zeros(nt)             # Mean elevation of ice-water interface
h_mean = np.zeros(nt)             # Mean elevation of surface above ice-water interface
x_left = np.zeros(nt)             # Left grounding line position
x_right = np.zeros(nt)            # Right grounding line position
P_res = np.zeros(nt)              # Penalty functional residual

if model == 'lake':
    lake_vol = np.zeros(nt)           # Lake volume
    dPw = np.zeros(nt)                # Water pressure - hydrostatic pressure

t = 0                             # time

# Begin time stepping
for i in range(nt):

    print('-----------------------------------------------')
    print('Iteration '+str(i+1)+' out of '+str(nt))

    if t==0:
        # Set initial conditions.
        s_mean_i = s_mean0                    # Mean ice-water elevation.
        F_h = lambda x: Hght                  # Ice-air surface function
        F_s = lambda x: interface(x)    # Lower surface function

        if model == 'marine':
            w = get_zero_m(mesh)              # Placeholder for first iteration.
        elif model == 'lake':
            w = get_zero_l(mesh)

        mesh,F_s,F_h,s_mean_i,h_mean_i,XL,XR = mesh_routine(w,mesh,dt)

    # Solve the Stoke problem.
    # Returns solutions "w" and penalty functional residual "Perr_i"
    if model == 'marine':
        w,P_res_i = stokes_solve_marine(mesh,F_h,t)

    elif model == 'lake':
        w,P_res_i = stokes_solve_lake(mesh,lake_vol_0,s_mean_i,F_h,t)

    # Solve the surface kinematic equations, move the mesh, and compute the
    # grounding line positions.
    mesh,F_s,F_h,s_mean_i,h_mean_i,XL,XR = mesh_routine(w,mesh,dt)

    # Save quantities of interest.
    P_res[i] = P_res_i
    s_mean[i] = s_mean_i
    h_mean[i] = h_mean_i
    x_left[i] = XL
    x_right[i] = XR
    Gamma_s[:,i] = F_s(X_fine)
    Gamma_h[:,i] = F_h(X_fine)

    # Save (u,p) solution for viewing in Paraview.
    if model == 'lake':
        # Save additional quantities unique to the lake problem.
        #Compute lake volume: Integral of lower surface minus the bed elevation    x_right[i] = XR.
        lake_vol[i] = scpint.quad(lambda x: F_s(x)-bed(x),0,Lngth,full_output=1)[0]

        #Compute difference between hydrostatic pressure and mean water pressure
        dPw[i] = (np.abs(w.sub(2).compute_vertex_values(mesh)[0])-rho_i*g*(h_mean_i-s_mean_i))/1.0e3

        # Save Stokes solution
        _u, _p,_pw = w.split()
        _u.rename("vel", "U")
        _p.rename("press","P")
        vtkfile_u << (_u,t)
        vtkfile_p << (_p,t)


    elif model == 'marine':
        # Save Stokes solution
        _u, _p = w.split()
        _u.rename("vel", "U")
        _p.rename("press","P")
        vtkfile_u << (_u,t)
        vtkfile_p << (_p,t)

    # Update time
    t += dt

    # Print information of interest.
    print('Left grounding line: '+str(x_left[i]/1000.0)+' km')
    print('Right grounding line: '+str(x_right[i]/1000.0)+' km')

    if model == 'lake':
        # Print information of interest unique to the lake problem.
        print('Water volume percent change: '+str(100*(lake_vol[i] - lake_vol[0])/lake_vol[0]))


# Save quantities of interest.
t_arr = np.linspace(0,t_final,num=int(nt_per_year*t_final/3.154e7))

np.savetxt('results/Gamma_s',Gamma_s)
np.savetxt('results/Gamma_h',Gamma_h)
np.savetxt('results/s_mean',s_mean)
np.savetxt('results/h_mean',h_mean)
np.savetxt('results/x_left',x_left)
np.savetxt('results/x_right',x_right)
np.savetxt('results/P_res',P_res)
np.savetxt('results/X',X_fine)           # X = spatial coordinate
np.savetxt('results/t',t_arr)            # t = time coordinate


if model == 'marine' and tides=='off':
    new_mesh << mesh

if model == 'lake':
    # Save quantities of interest unique to the lake problem.
    np.savetxt('results/dPw',dPw)
    np.savetxt('results/lake_vol',lake_vol)



#----------------------------PLOTTING-------------------------------------------

# For the test problem and paper problems only

#------------------------- Plotting for test problem----------------------------
if model_setup == 'test' :
    testplot_1(s_mean,h_mean,lake_vol,x_right,x_left,dPw)
    testplot_2(Gamma_s,Gamma_h,x_left,x_right)

#------------------ Plotting for Figures 2-3 in paper---------------------------
if model_setup == 'tides_paper':
    paperplot_fig2(s_mean,h_mean,x_right,x_left)
    paperplot_fig3(Gamma_s,Gamma_h,x_left,x_right)


#------------------ Plotting for Figures 4-5 in paper---------------------------
# Also plots the graphical abstract (gabstract.png)
if model_setup == 'lake_paper':
    paperplot_fig4(s_mean,h_mean,lake_vol,x_right,x_left,dPw)
    paperplot_fig5(Gamma_s,Gamma_h,x_left,x_right)
    gabstract(Gamma_s,Gamma_h,x_left,x_right)
