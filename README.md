# 2.5 Layer Reduced Gravity Model
This code is for a 2.5 layer linearized reduced gravity model - i.e. 2 active layers over a quiesscent abyss. It is forced by a zonal wind only through a channel that the user specifies. The boundaries are periodic (re-entrant channel) where the zonal wind is located and otherwise there is no normal flow through boundaries and a no slip condition applied. The model is built on an Arakawa C grid.

Much of the model code structure is a modified version of the single layer linearized shallow water code by James Penn (see https://empslocal.ex.ac.uk/people/staff/gv219/codes/linearshallowwater.py)

This repo is divided such that the grid functions and bulk of the model are located in model python files (e.g., model_force_north_sponge_rk4_noh_diff_noslip_properimplement_update_onlycorners_residualcirculation_nonlinearcontinuity.py). The model can be run in execute python files found in each subfolder which corresponds to a parameter set-up.

This main model file has been updated and has the following features (as of March 2024):

    Nonlinear continuity equations, linear momentum equations
    Can add mass flux at the southern boundary, which can be given as a single value (implemented as a step function) or as a vector (if time varying)
    Optional mass transfer between layers can be prescribed
    Sponge layers which damp just velocity (not h)
    On beta plane
    RK4 solver in time
    No slip and no normal flow boundary conditions

The most essential part of the directory is in fixed_boundary_rk4_novolumeleak/finalized_match_mitgcmdomainsize/ which is the 2.5 layer reduced gravity model in its most recent form, with domain size set to that of MITgcm. We also test our results against the full 3 layer model, which is in fixed_boundary_rk4_novolumeleak/3layer_check

## Figures
The figures in the paper are made in the following files. 

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

2. The subfolder full_cor_version has a full coriolis expression rather than beta plane (although still Cartesian grid), but this is a much older version of the work which was using free-slip boundary conditions and an euler time step and is <strong>NOT</strong> used for any results in the manuscript. It is kept here so that one can see how to modify the more recent model file to use a full Coriolis expression rather than beta plane.
