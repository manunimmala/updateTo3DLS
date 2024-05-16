#!/bin/bash

#SBATCH --job-name LS
#SBATCH --account=esm_5554
#SBATCH --partition=normal_q
#SBATCH --cpus-per-task=96
#SBATCH -t 16:00:00
#SBATCH --output=LS.out
#SBATCH --error=LS.err

#SBATCH --mail-type=all          
#SBATCH --mail-user=nimmala@vt.edu

module reset
module load tinkercliffs-rome/matlab/R2021a
matlab -batch launch_LS

exit 0
