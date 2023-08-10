#!/usr/bin/bash --login
#SBATCH --job-name=Permu7
#SBATCH --time=8:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-125
#SBATCH --constraint=[intel16|intel18]
#SBATCH --mail-user=smit2276@msu.edu
#SBATCH --mail-type=ALL

cd $SLURM_SUBMIT_DIR

R CMD BATCH --no-save --no-restore '--args trait="bhmd" genotypes="calls" jobID='${SLURM_ARRAY_TASK_ID} 1_Permu_7.R  out_${SLURM_ARRAY_TASK_ID}

# This is a 'slurm' submission sctipt, I've included this because otherwise 
# 1_Permu_7.R makes less sense
#  
# jobID is defined here, 
# where ${SLURM_ARRAY_TASK_ID} is 1 through 125 for each job. 
#
# out_* file is where the teminal output goes, including error messages
# if there is an error like going over the time limit, 
# or if the job is finished, my MSU email will recive an alert.
# 
