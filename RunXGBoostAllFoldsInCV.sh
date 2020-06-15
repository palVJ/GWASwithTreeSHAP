#!/bin/sh
#SBATCH --mem=120G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-4
#SBATCH --cpus-per-task=20
#SBATCH --time=03:00:00
#SBATCH --mail-user=youremail

module purge
module load miniconda
source activate r_env


/usr/bin/time -v Rscript --vanilla RunXGBoostInCVOnSubset.R ${SLURM_ARRAY_TASK_ID} ${shID} ${eta} ${colsample_bylevel} ${subsample} ${colsample_bytree} ${max_depth} 
${nthreads}




