#!/bin/bash

#SBATCH --job-name=tar
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10GB
#SBATCH --time=0:40:00
#SBATCH --mail-type=END
#SBATCH --mail-user=abf376@nyu.edu
#SBATCH --output=slurm_%j.out

singularity exec --overlay /scratch/$USER/my_env/overlay-5GB-200K.ext3:ro /scratch/work/public/singularity/centos-8.2.2004.sif /bin/bash -c "source /ext3/env.sh; python -u /scratch/abf376/2.5_layer_model/3layer_check/execute_files/execute_layer1_masspert_fromstationary_rk4_visc8E3_noslip_equator0_realthreelayer_lowerlayer10000_cuttimestep.py"