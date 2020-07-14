# All model parameters and options are recorded here.
import numpy as np
#-------------------------------------------------------------------------------
#-----------------------------MODEL OPTIONS-------------------------------------

# Default parameters correspond to the test case, which outpouts png plots
# that should be compared to the plots in the 'test results' folder

model_setup = 'test'                        # set for quick example

#model_setup = 'wedge_test'                 # set for 'wedge' tests

# Set model to 'lake' for subglacial lake or 'marine' for marine ice sheet:
model = 'lake'

# For the marine ice sheet setup, turn the tidal cycle 'on' or 'off'.
# (Turn off for lake problem)

tides = 'off'

# Turn 'on' or 'off' real-time plotting that saves a png figure called 'surfs' at
# each time step of the free surface geometry.

realtime_plot = 'on'

# Turn 'on' or 'off' Newton convergence information:
print_convergence = 'off'

# Mesh resolution at the lower boundary
DX_s = 200.0                    # Element width at lower boundary (in meters)
                              # Default values are {200,100,50,25,12.5}
                              # This is used for (1) setting the element width in
                              # gendomain.py and (2) selecting the mesh in main.py.

DX_h = 250.0                  # Element width at the upper surface (in meters)


if model != 'marine' and model != 'lake':
    sys.exit('ERROR: Set \'model\' to \'marine\' or \'lake\' ONLY in params.py file!')

#-------------------------------------------------------------------------------


#-----------------------------MODEL PARAMETERS----------------------------------
#-------------------------------------------------------------------------------

# Material parameters
A0 = 3.1689e-24                    # Glen's law coefficient (ice softness)
n = 3.0                            # Glen's law exponent
rm2 = 1 + 1.0/n - 2.0              # Exponent in variational forms: r-2
B0 = A0**(-1/n)                    # Ice hardness
B = (2**((n-1.0)/(2*n)))*B0        # "2*Viscosity" constant in weak form
rho_i = 917.0                      # Density of ice
rho_w = 1000.0                     # Density of water
g = 9.81                           # Gravitational acceleration
C = 1.0e5                          # Sliding law friction coefficient

# Numerical parameters
eps_v = 1.0e-15                    # Flow law regularization parameter
eps_p = 1.0e-13                    # Penalty method parameter
quad_degree = 16                   # Quadrature degree for weak forms

tol = 1.0e-3                       # Numerical tolerance for boundary geometry:
                                   # s(x,t) - b(x) > tol on ice-water boundary,
                                   # s(x,t) - b(x) <= tol on ice-bed boundary.

# Geometry parameters
Lngth = 10*1000.0                  # Length of the domain
Hght = 1000.0                      # (Initial) Height of the domain

sea_level = Hght*(917.0/1000.0)    # Sea level elevation.
                                   # (Initial sea level for the tides problem)
# Time-stepping parameters for non-tidal problems
nt_per_year = 250                  # Number of timesteps per year.
t_final =1*3.154e7                 # Final time (yr*sec_per_year).

nt = int(nt_per_year*t_final/3.154e7) # Number of time steps
dt = t_final/nt                       # Timestep size

nx = 10*Nx
X_fine = np.linspace(0,Lngth,num=nx)  # Horizontal coordinate for computing surface
                                      # slopes and plotting.

# Set zero inflow/outflow speed boundary condition for subglacial lake problem
U0  = 0.0                             # Zero inflow/outflow speed

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
