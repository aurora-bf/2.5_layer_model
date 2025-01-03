# Code for ''A theory for how the depth of meltwater injection impacts regional sea level evolution"
Preprint: https://essopenarchive.org/users/751256/articles/1215851-an-analytical-theory-for-the-sensitivity-of-regional-sea-level-adjustment-to-the-depth-of-antarctic-meltwater-fluxes

This repo has code for the numerical 2.5-layer model utilized throughout the paper, implementation of the analytic solution, and figure creation.
## Numerical 2.5 layer reduced gravity model
The numerical model is a 2.5 layer linearized reduced gravity model - i.e. 2 active layers over a quiesscent abyss. It is forced by a zonal wind only through a channel that the user specifies. The boundaries are periodic (re-entrant channel) where the zonal wind is located and otherwise there is no normal flow through boundaries and a no slip condition applied. The model is built on an Arakawa C grid.

Much of the model code structure is a modified version of the single layer linearized shallow water code by James Penn (see https://empslocal.ex.ac.uk/people/staff/gv219/codes/linearshallowwater.py)

The numerical experiments are conducted such that the grid functions and bulk of the model are located in model python files (e.g., model_force_north_sponge_rk4_noh_diff_noslip_properimplement_update_onlycorners_residualcirculation_nonlinearcontinuity.py). The model can be run in execute python files found in each subfolder which corresponds to a parameter set-up.

The main model file has the following features:

    Nonlinear continuity equations, linear momentum equations
    Can add mass flux at the southern boundary, which can be given as a single value (implemented as a step function) or as a vector (if time varying)
    Optional mass transfer between layers can be prescribed
    Sponge layers which damp just velocity (not h)
    On beta plane
    RK4 solver in time
    No slip and no normal flow boundary conditions


The key part of this directory is /main_2.5layer/. Here, you will find the model file, execute files as necessary in subfolders corresponding to different parameter set-ups and notebooks analyzing the resutls (see description of the notebooks in the Figures section below).

## Figures
The figures in the paper are made in notebooks as listed below. In these notebooks, as necessary, the analytic solution is computed, compared to the numerical runs, and plots are created. 

fixed_boundary_rk4_novolumeleak/finalized_match_mitgcmdomainsize/sealevel_fixeddomainsize_250m_rhodiff2.ipynb: Figs 3, 4, 5, 6, A1, B1

fixed_boundary_rk4_novolumeleak/finalized_match_mitgcmdomainsize/domain_diagram.ipynb: Fig 2

fixed_boundary_rk4_novolumeleak/finalized_match_mitgcmdomainsize/sealevel_fixeddomainsize_varyparameters.ipynb: Fig 7

fixed_boundary_rk4_novolumeleak/finalized_match_mitgcmdomainsize/sealevel_fixeddomainsize_250m_rhodiff2_unequalstratification.ipynb: Figs 8, B2

fixed_boundary_rk4_novolumeleak/finalized_match_mitgcmdomainsize/sealevel_fixeddomainsize_250m_rhodiff2_wind_morespunup.ipynb: Fig 9

fixed_boundary_rk4_novolumeleak/3layer_check/results_3layer_cleaned_updatedparameter.ipynb: Suppplementary figure


Fig 1 is a schematic not made in this repo.


As an example (not used much in the paper other than as part of Fig 7), some analogues for a different parameter set-up are included. For example, equivalents to Figs 1,4,5,6,D1,F1 are made in sealevel_fixeddomainsize_100m_rhodiff2.ipynb. 

## Additional notes 
1. The conda environment used when creating this repo is in environment_apr2024.yml and the environment can be recreated using conda env create -f environment_apr2024.yml

2. The 3layer_check subfolder is not used to make any figures in the manuscript. It is merely to show a sanity check, that our results are true with a three layer model (with deep abyss), in addition to a 2.5-layer model, as we have focused on.
