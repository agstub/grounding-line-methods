[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4302610.svg)](https://doi.org/10.5281/zenodo.4302610)

grounding-line-methods

Author: Aaron Stubblefield (Columbia University, LDEO).

# Overview
This repository contains FEniCS python code for simulating grounding line migration in
the subglacial lake and marine ice sheet settings. The model is
isothermal Stokes flow with nonlinear ("Glen's law") viscosity. The contact
conditions that determine whether ice remains in contact with the bed or goes
afloat are enforced with a penalty functional. A full description of the model
formulation is described in a manuscript submitted to JFM.

# Dependencies
## Required dependencies
As of this commit, this code runs with the latest FEniCS Docker image (https://fenicsproject.org/download/).
Docker may be obtained at: https://www.docker.com/. To run the Docker image:

`sudo docker run -ti -p 127.0.0.1:8000:8000 -v $(pwd):/home/fenics/shared -w /home/fenics/shared quay.io/fenicsproject/stable:latest`


## Optional dependencies

1. Gmsh (http://gmsh.info/) is used for mesh generation. Note that the default meshes
used in the manuscript are already included in the repository (*meshes* directory, see below).

2. FFmpeg (https://www.ffmpeg.org/) can be used, along with **make_movie.py**,
to create a video of the evolving free surface geometry over time. See description below.

3. ParaView is useful for visualizing velocity/pressure solutions to the Stokes equations (https://www.paraview.org/).

# Contents

## 1. Source files
The model is organized in 7 python files in the *source* directory as follows.

1. **geometry.py** contains the geometric description of the bed and initial ice-water interface.

2. **params.py** contains all of the model parameters and model options.

3. **stokes.py** contains the Stokes system solver and related functions.

4. **meshfcns.py** contains functions that solve the surface kinematic equations, move the mesh,
    and return the grounding line positions.

5. **boundaryconds.py** contains functions that mark the mesh boundary and apply boundary conditions.

6. **hydrology.py** contains the subglacial lake volume change
timeseries and the sea level change timeseries.

7. **main.py** runs the model. It contains the time-stepping loop that
calls the Stokes solver and mesh-related functions at each timestep, and saves the output.

The *model_examples* subdirectory contains alternative parameter files
that produce the results in the manuscript. See the 'Running problems from the paper'
section below.

## 2. Scripts

The *scripts* directory contains:

1. **create_mesh.py**: generates meshes using **geometry.py** and Gmsh.

2. **plotting.py**: When the default options are used, this is called in **main.py** to creates the 'test' plots (laketest1.png and laketest2.png). These should match the .png's in the *testresults* directory. When
**params.py** is replaced by either **params_tides.py** or **params_lake.py**
from *source/model_examples*, this function is called in **main.py** and the relevant figures from the paper are produced.

3. **make_movie.py**: generates .png
images of the basal and upper surfaces for each time step. These may then be
used to create a .mp4 movie using the FFmpeg command in
the comments at the top of the **make_movie.py** file:
`ffmpeg -r frame_rate -f image2 -s 1920x1080 -i %01d.png -vcodec libx264 -pix_fmt yuv420p -vf scale=1280:-2 movie.mp4`
where `frame_rate` is an integer (e.g., 50).

## 3. Meshes
The *meshes* directory contains the .xml mesh files used for the simulations in the paper.
All meshes have an element
width of 250 m at the upper surface. The element width at the lower surface
is different for each (model type)_DXp.xml file, where the integer p is the element width at the
lower surface in meters (rounded down to 12 in the case DX=12.5).
The tides_DX12.xml and marine_DX12.xml meshes have element widths of 12.5 m at the lower
surface.

## 4. Test results
The *test_results* directory contains two .png images (lake_test_1.png and lake_test_2.png) that are produced when
the model is run with the default parameters.

# Running the test problem
The test problem simulates a subglacial lake filling-draining cycle with a period
of one year. The mesh spacing at the lower boundary is 200 m and the timestep size
is 1/250 yr. These are the default parameters in **params.py**. These parameters
are chosen because the run time is short--about 4.5 minutes on a Dell XPS-13 (9360).

To run the test problem:

1. Run the FEniCS Docker image.

2. In Docker, run the main file from the parent directory: `python3 ./source/main.py`

Upon completion, the code should produce two .png files (lake_test_1.png and lake_test_2.png) that match the .png's in the *test_results* directory.

# Running problems from the paper
To reproduce Figures 2-3 in the paper, replace **params.py** with
**params_tides.py** from *source/model_examples* (renamed to **params.py**) and run the command
in the previous section.

To reproduce Figures 5-6 in the paper, replace **params.py** with
**params_lake.py** from *source/model_examples* (renamed to **params.py**) and run the command
in the previous section.

The convergence results can be obtained by modifying
`DX_s`, `eps_p`, and `nt_per_year` in the **params.py** file as specified in the
appendix of the manuscript. The "wedge test" examples can be run by setting
`model_setup="wedge_test"`. See **params_wedge.py** in the *source/model_examples* directory.

# Output

Model output is saved in a *results* directory. This includes

1. the Stokes solution vtk files (*stokes* subdirectory),

2. upper and lower surface geometry at each timestep (Gamma_s and Gamma_h, resp.),

3. minimum and maximum grounding line positions (x_left and x_right, resp.),

4. penalty functional residual (P_res),

5. mean elevation of upper and lower surfaces over the ice-water interface (s_mean and h_mean, resp.),

6. spatial coordinate (X), and

7. time coordinate (t).

    For the subglacial lake model, the

8. lake volume timeseries (lake_vol) and

9. deviation of mean water pressure from cryostatic pressure (dPw, in kilopascals)

are also saved.
Gamma_s and Gamma_h are two-dimensional arrays:
the columns are the free surfaces at the timestep corresponding to the row index.
The spatial coordinate X is finer than the mesh spacing because Gamma_s and Gamma_h
are created by linear interpolation of the mesh nodes in SciPy.

# Model options

Model options are set in the **params.py** file. For example:

1. Choose the *subglacial lake* or *marine ice sheet* model setups by setting
`model='lake'` or `model='marine'` at the top of the **params.py** file.

2. For the marine ice sheet setup, ocean tides may be turned on by setting
`tides='on'`.

3. *Real-time plotting* is available by setting `realtime_plot = "on"` in the
**params.py** file. This outputs a png called 'surfaces' of the free surface geometry at each
timestep.

Other parameters of interest in **params.py** that may be modified are `t_final` (final time)
`nt_per_yer` (number of timesteps per year), `eps_p` (penalty method parameter),
and `U0` (inflow/outflow speed).

The sea level change and lake volume change timeseries may be modified
in the **hydrology.py** file.

# Domain description and mesh generation
The bed geometry can be changed by modifying the `bed` function in **geometry.py**. This requires defining
a consistent initial lower boundary by modifying the `interface` function. A
new mesh must then be generated by running `python ./scripts/create_mesh.py` (requires Gmsh).
The element widths at the lower and upper boundaries are set by `DX_s` and `DX_h` in **params.py**, respectively.
The `mesh_name.msh` file generate by Gmsh is then converted to a .xml file via
`dolfin-convert mesh_name.msh mesh_name.xml` in the Docker image.

To generate the mesh for the tidal problem (i.e., tides_DX12 in *meshes* directory)
run the model with the parameter file **params_marine.py** (stored in *source/model_examples*).
Upon completion, this will create a new mesh called 'tides_DX(spacing).xml" that provides
a good initial geometry for tidal simulations.

# Acknowledgements
Many thanks to Ed Bueler (UAF) who introduced me to free surface Stokes problems during UAF's 2018 International Summer School in Glaciology, and to the UAF Glaciers group (https://glaciers.gi.alaska.edu/) for hosting the summer school. Marc Spiegelman and Tim Creyts (Columbia University, LDEO) advised me throughout this project.
