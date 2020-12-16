#!/bin/sh
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --ntasks=1
#SBATCH --array=1-470
#SBATCH --nodes=1
#SBATCH --time=00:22:00
#SBATCH --mail-user=ALL

module load miniconda
source activate r_env

/usr/bin/time -v Rscript --vanilla ~/ShapInteractions.R ${SLURM_ARRAY_TASK_ID}
