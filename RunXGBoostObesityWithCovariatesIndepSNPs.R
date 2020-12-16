#Run XGBoost in a cross-validation setting on a subset of all data consisting of SNPs and environmental data.
#The SNPs should have low mutual correlation.
#Run with the following arguments:
#The name o

library(roxygen2)
library(devtools)
document()
setwd("~/XGBoost/XGBoostOnSNPData")

args <- commandArgs(trailingOnly = TRUE)
WhichIsTestdata = as.numeric(args[1])

shID = args[2]

#Include covariates to obesity data:
library(data.table)
CovData = fread("~/ObesityCovOnlyCaucasiansIndepQC.txt",header = TRUE, sep  = "\t", colClasses = c("IID" = "character"))
#Which IIDs to include.
IIDsToInclude = fread(paste("~/RawFiles/IIDsToIncludeIndepSNPs",shID,".txt",sep = ""),header = TRUE, sep = "\t",
select = c("IID"),colClasses = c("IID"="character"))
CovData = subset(CovData,CovData$IID %in% as.character(IIDsToInclude$IID))

RUN = XGBoostOnSNPdata_cv(TrainData = c(paste("~/RawFiles/SNP_add_","SubsetIndepSNPs",shID,"_Trai1.raw",sep = ""),
paste("/SAY/dbgapstg/dbgap/work/paalvj/RawFiles/SNP_add_","SubsetIndepSNPs",shID,"_Trai2.raw",sep = ""),
paste("/SAY/dbgapstg/dbgap/work/paalvj/RawFiles/SNP_add_","SubsetIndepSNPs",shID,"_Trai3.raw",sep = ""),
paste("/SAY/dbgapstg/dbgap/work/paalvj/RawFiles/SNP_add_","SubsetIndepSNPs",shID,"_Trai4.raw",sep = "")),
ValData = paste("/SAY/dbgapstg/dbgap/work/paalvj/RawFiles/SNP_add_","SubsetIndepSNPs",shID,"_ValData.raw",sep = ""),
Kround = WhichIsTestdata,hyperparameters = list(eta = 0.05,max.depth = 2,colsample_bytree = 0.8,colsample_bylevel=0.8,subsample = 0.8,tree_method = "hist",
grow_policy = "lossguide"),nrounds = 8000,early_stopping_rounds = 20,nthreads = 20,covariatesTable = CovData)


#Write importance table to file:
fwrite(RUN$Importance,paste("~/ImportanceTables/BMI_ImportanceUnifWholeGenomeWithCovsIndepSNPs","Subset",shID,"TestData",WhichIsTestdata,".txt",
sep = ""),quote = FALSE, row.names = FALSE, sep = "\t",nThread=2)


#Save R-object models:
save(RUN, file = paste("~/XGBmodels/XGB_UnifWholeGenomeWithCovsIndepSNPs","Subset",shID,"Testdata",WhichIsTestdata,".RData",sep = ""))
