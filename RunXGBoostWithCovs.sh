#!/bin/sh
#SBATCH --mem=140G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-4
#SBATCH --cpus-per-task=22
#SBATCH --time=10:00:00
#SBATCH --mail-user=pal.johnsen@yale.edu

module purge
module load miniconda
source activate r_env


/usr/bin/time -v Rscript --vanilla RunXGBoostObesityWithCovariates.R ${SLURM_ARRAY_TASK_ID} "$shID"
