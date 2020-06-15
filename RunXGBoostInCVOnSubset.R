#Run XGBoost based on a sample subset of individuals with uncorrelated SNPs used in the ranking process.

library(roxygen2)
library(devtools)
#Put folder destination of XGBoostOnSNPData here:
setwd("~/XGBoostOnSNPData")
document()

args <- commandArgs(trailingOnly = TRUE)
#WhichIsValidData is the fold in cross-validation that is used to validate during training.
WhichIsValidData = as.numeric(args[1])
shID = args[2]
eta = as.numeric(args[3])
colsample_bylevel = as.numeric(args[4])
subsample = as.numeric(args[5])
colsample_bytree = as.numeric(args[6])
max_depth = as.numeric(args[7])
nthreads = as.numeric(args[8])

#TrainData is the character vector of filenames for all folds, including validation data, that is used for training in cross-validation.
#TestData is the fold used to test the model created from the TrainingData.


RUN = XGBoostOnSNPdata_cv(TrainData = c(paste("~/SNP_add_Subset",shID,"_Trai1.raw",sep = 
""),paste("~/SNP_add_Subset",shID,"_Trai2.raw",sep = ""),
paste("~/SNP_add_Subset",shID,"_Trai3.raw",sep = ""),paste("~/SNP_add_Subset",shID,"_Trai4.raw",sep = "")),
TestData = paste("~/SNP_add_Subset",shID,"_TestData",".raw",sep = ""),
Kround = WhichIsValidData,hyperparameters = list(eta = eta,max_depth = max_depth,colsample_bytree = colsample_bytree,subsample = subsample,colsample_bylevel = 
colsample_bylevel,tree_method = "hist",grow_policy= "lossguide"),
nrounds = 8000,early_stopping_rounds = 20 ,nthreads = nthreads)

#Save R-object models:
save(RUN, file = paste("~/XGB_UnifWholeGenomeUnCorSNPs","Subset",shID,"ValidData",
WhichIsValidData,"e",eta,"ctscl",colsample_bytree,subsample,colsample_bylevel,"md",max_depth,".RData",sep 
= ""))
