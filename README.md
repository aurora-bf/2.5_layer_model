# Code for ''A theory for how the depth of meltwater injection impacts regional sea level evolution"
Paper: https://journals.ametsoc.org/view/journals/phoc/aop/JPO-D-24-0153.1/JPO-D-24-0153.1.xml

This repo has code for the numerical 2.5-layer model utilized throughout the paper, implementation of the analytic solution, and figure creation.
## Numerical 2.5 layer reduced gravity model
The numerical model is a 2.5-layer linearized reduced gravity model - i.e. 2 active layers over a quiesscent abyss. It is (optionally) forced by a zonal wind through a channel that the user specifies. The boundaries are periodic in the re-entrant channel (where the zonal wind is located) and otherwise there is no normal flow through boundaries and a no slip condition applied. The model is built on an Arakawa C grid.

**Much of the model code structure is a modified version of the single layer linearized shallow water code by Penn and Vallis** (see https://empslocal.ex.ac.uk/people/staff/gv219/codes/linearshallowwater.py). The equations solved are different and we have added different boundary conditions + an RK4 solver as well as other modifications. However, some functions/comments remain unmodified as they originally appeared in the cited code above.

The numerical experiments are conducted such that the grid functions and bulk of the model are located in model python files (e.g., /main_2.5layer/model.py). The model can be run in execute python files found in each subfolder which corresponds to a parameter set-up (e.g., /main_2.5layer/250m_rhodiff2/execute_layer1_masspert_fromstationary_rk4_visc8E3_noslip_equator0_nonlinear_250m_g1equalg2.py).

The main model file has the following features:

    Nonlinear continuity equations, linear momentum equations
    Can add mass flux at the southern boundary, which can be given as a single value (implemented as a step function) or as a vector (if time varying)
    Optional mass transfer between layers can be prescribed
    Sponge layers which damp velocity (not height)
    On beta plane
    RK4 solver in time
    No slip and no normal flow boundary conditions, except in re-entrant channel


The key part of this directory is /main_2.5layer/. Here, you will find the model file, execute files as necessary in subfolders corresponding to different parameter set-ups, and notebooks analyzing the results (see description of the notebooks in the Figures section below).

## Figures
The figures in the paper are made in notebooks as listed below. In these notebooks, as necessary, the analytic solution is computed, compared to the numerical runs, and plots are created. 

main_2.5layer/analytic_defaultparameters.ipynb: Figs 3, 4, 5, A1, B1, B2

main_2.5layer/domain_diagram.ipynb: Fig 2

main_2.5layer/varyparameters.ipynb: Fig 6

main_2.5layer/plotting_figure7.ipynb: Fig 7. This main notebook depends on results produced in analytic_unequalstratification.ipynb and analytic_wind.ipynb.
  

Fig 1 is a schematic not made in this repo.

Using numerical results for different parameter set-ups (used in Figure 6), one could easily modify the above notebooks to get analytic results for different set-ups as well.

Large files produced in execute files and needed for recreating the analytic results are available at: https://drive.google.com/drive/folders/1Z1Tk683nuRAuMG4KPhCP8TACMFuEWxJC?usp=sharing so one doesn't need to run the supplied execute files to recreate the main results in the paper.

## Additional notes 
1. The conda environment used when creating this repo is in environment_apr2024.yml and the environment can be recreated using conda env create -f environment_apr2024.yml

2. The 3layer_check subfolder is not used to make any figures in the manuscript. It is merely to show a check that we tested that our numerical results found in the 2.5-layer model are also true with a 3-layer model with a deep abyss.

