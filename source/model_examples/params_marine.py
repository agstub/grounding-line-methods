# All model parameters and options are recorded here.

# NOTE: These parameters correspond to the marine problem in the paper.
# This is used to construct the initial geometry for the tidal problem.

import numpy as np
#-------------------------------------------------------------------------------
#-----------------------------MODEL OPTIONS-------------------------------------

# flag for constructing mesh for the tidal problem
# (saves mesh at final time of simulation)
model_setup = 'marine_paper'

# Set model to 'lake' for subglacial lake or 'marine' for marine ice sheet:

model = 'marine'

# For the marine ice sheet setup, turn the tidal cycle 'on' or 'off'.
# Turning tides 'on' defaults to a simulation time of ~1 week:

tides = 'off'

# Turn 'on' or 'off' real-time plotting that saves a png figure called 'surfs' at
# each time step of the free surface geometry.

realtime_plot = 'off'

# Turn 'on' or 'off' Newton convergence information:
print_convergence = 'on'

# Mesh resolution at the lower boundary
DX_s = 12.5                   # Element width at lower boundary (in meters)
                              # Default values are {200,100,50,25,12.5}
                              # This is used for (1) setting the element width in
                              # gendomain.py and (2) selecting the mesh in main.py.

DX_h = 250.0                  # Element width at the upper surface (in meters)


if model != 'marine' and model != 'lake':
    sys.exit('ERROR: Set \'model\' to \'marine\' or \'lake\' ONLY in params.py file!')

#-------------------------------------------------------------------------------


#-----------------------------MODEL PARAMETERS----------------------------------
#-------------------------------------------------------------------------------
# physical units:
# time - seconds
# space - meters
# pressure - pascals
# mass - kg

# Material parameters
A0 = 3.1689e-24                    # Glen's law coefficient (ice softness, Pa^{-n}/s)
n = 3.0                            # Glen's law exponent
rm2 = 1 + 1.0/n - 2.0              # Exponent in variational forms: r-2
B0 = A0**(-1/n)                    # Ice hardness (Pa s^{1/n})
B = (2**((n-1.0)/(2*n)))*B0        # "2*Viscosity" constant in weak form (Pa s^{1/n})
rho_i = 917.0                      # Density of ice (kg/m^3)
rho_w = 1000.0                     # Density of water (kg/m^3)
g = 9.81                           # Gravitational acceleration (m/s^2)
C = 1.0e5                          # Sliding law friction coefficient (Pa s^{1/n}/m)

# Numerical parameters
eps_v = 1.0e-15                    # Flow law regularization parameter
eps_p = 1.0e-13                    # Penalty method parameter
quad_degree = 16                   # Quadrature degree for weak forms

tol = 1.0e-3                       # Numerical tolerance for boundary geometry:
                                   # s(x,t) - b(x) > tol on ice-water boundary,
                                   # s(x,t) - b(x) <= tol on ice-bed boundary.

# Geometry parameters
Lngth = 20*1000.0                  # Length of the domain (m)
Hght = 500.0                       # (Initial) Height of the domain (m)

sea_level = Hght*(917.0/1000.0)    # Sea level elevation (m).
                                   # (Initial sea level for the tides problem)

# Time-stepping parameters for tidal problems
nt_per_year = 2000                 # Number of timesteps per year.
t_final = 0.35*3.154e7             # Final time (yr*sec_per_year).

nt = int(nt_per_year*t_final/3.154e7) # Number of time steps
dt = t_final/nt                       # Timestep size

nx = 10*int(Lngth/DX_s)               # Horizontal coordinate for computing surface
X_fine = np.linspace(0,Lngth,num=nx)  # slopes, interpolated grounding line positions, and plotting.

# Set positive inflow speed boundary condition for marine ice sheet problem
U0   = 1000.0/3.154e7                 # Inflow speed 1 km/yr (m/yr / sec/yr)

#-------------------------------------------------------------------------------
