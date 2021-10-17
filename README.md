# GWASwithTreeSHAP
A new method to explore gene-gene and gene-environment interactions in a GWAS using SHAP values.

## Published paper
The paper is published in BMC Bioinformatics, and is available here: https://doi.org/10.1186/s12859-021-04041-7

## Tutorial
In this section, we will present a tutorial for how to apply the code given GWAS data.
The code is constructed for use in a high performance computing (HPC) system using slurm.
The code is written in R as well as plink.
R-packages required are data.table and xgboost.

Before applying the code, first separate your GWAS data in three disjoint data sets:
The ranking data, the fitting data and the evaluation data.

### Ranking process
In the ranking process we apply the ranking data, of size N, for ranking the features by importance, and thereby removing noise in the later process.
We create A randomly chosen subsets consisting of G <= N individuals and S SNPs with small mutual correlation. This can be done by running
the batch script **MakeCrossValdsIndepSNPs.sh** in the XGBoostAndTreeSHAP folder, with each line producing one subset. Add or remove lines for your chosen
number of A subsets. (Alternatively, make a batch script of A array jobs running **PrepareGWASdataIndepSNPs.sh**). Before running this file, you will need to
include a file of IIDs of the individuals in the ranking data, as well as a file including all SNPs available in **SampleSNPAndIndsUniformly.R**, as well as parameters
explained in detail below.

In the batch script **PrepareGWASdataIndepSNPs.sh**, you include the name of the genotype file (including .bed, .bim and .fam files).
Here, you also decide the parameters G, as well as a parameter H, using the R-function **SampleSNPAndIndsUniformly.R**, indicating the number of SNPs to apply in the function
[--indep-pairwise](https://www.cog-genomics.org/plink/1.9/ld) in plink. You will need to tune this parameter H, in order to eventually get around S SNPs with small mutual correlation.
Here, you also decide the parameters in --indep-pairwise such as r^2 (pearson correlation squared), the window size and the step size.
Restricting the minimal MAF of the SNPs as well as the minimum p-value from the Hardy-Weinberg equilibrium exact test p-value, can be set in the plink code provided in
**PrepareGWASdataIndepSNPs.sh**.

The function **SplitDataForCV.R** further splits the data set in K sets to prepare for cross-validation. The final data sets
are [.raw](https://www.cog-genomics.org/plink/1.9/formats#raw) files, where the SNPs by default are in additive form (allele counts from 0 to 2).
If you want the allele counts in another form, change the parameter --recode A to one of the others [options](https://www.cog-genomics.org/plink/1.9/data#recode).
Save all data files in one directory (such as RawFiles). Let the saved files be in the form "SNP_add_SubsetIndepSNPs**shID**_Trai**i**.raw", with **shID** indicating an ID
of the in total A subsets, and **i**, the ID one of the folds used during cross-validation.

Run the function **IndepCrossValdsInRankingProcess.sh** which includes the batch script **RunXGBoostWithCovs.sh** which construct array jobs of size K-1 that runs
the R-function **RunXGBoostObesityWithCovariates.R**. This R-function requires the R-packages xgboost and data.table. For instance let these packages be downloaded in your own conda R environment
r_env. In **RunXGBoostObesityWithCovariates.R**, include the environmental covariates to use in CovData. The hyperparameters of the constructed xgboost-models is set here as an argument
in the function **XGBoostOnSNPdata_cv**.
Create a directory consisting of all the constructed xgboost-models, and save them in the form "XGB_UnifWholeGenomeWithCovsIndepSNPsSubset**shID**Testdata**WhichIsTestdata**.RData",
where **WhichIsTestdata** indicates which folder that is not used during the construction of the model in the cross-validation. 
