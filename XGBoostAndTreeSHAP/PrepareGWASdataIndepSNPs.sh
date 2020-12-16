#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G
#SBATCH --nodes=1
#SBATCH --time=00:50:00
#SBATCH --mail-user=ALL

module purge
module load R

Rscript ~/SampleSNPAndIndsUniformly.R 157000 70000 "$shID"


module load PLINK/1.90-beta5.3

plink --bfile ~/genotype_files/pleiotropy_geneticfiles/UKB_Caucasians_phenotypeindepqc120319 \
--extract ~/RawFiles/snpsToIncludeBeforeIndep${shID}.txt \
--keep ~/RawFiles/IIDsToIncludeIndepSNPs${shID}.txt \
--indep-pairwise 50 5 0.2 --out ~/RawFiles/IndepSNPsToInclude${shID}

plink --bfile ~/genotype_files/pleiotropy_geneticfiles/UKB_Caucasians_phenotypeindepqc120319 \
--pheno ~/PhenoFileOnlyCaucasiansIndsQC.txt \
--1 \
--pheno-name "Obese" \
--recode A \
--maf 0.01 \
--hwe 5e-8 \
--extract ~/RawFiles/IndepSNPsToInclude${shID}.prune.in \
--keep ~/RawFiles/IIDsToIncludeIndepSNPs${shID}.txt \
--out ~/RawFiles/SNP_additive_data_unifIndepSNPs${shID} \

Rscript ~/SplitDataForCV.R 5 ${shID}
rm ~/RawFiles/snpsToIncludeBeforeIndep${shID}.txt
rm ~/RawFiles/IIDsToIncludeIndepSNPs${shID}.prune.out
rm ~/RawFiles/IIDsToIncludeIndepSNPs${shID}.log
rm ~/RawFiles/SNP_additive_data_unifIndepSNPs${shID}.raw
