#Compute relative contribution of interactions for 100 people

#Array-ID
args = commandArgs(trailingOnly=TRUE)
arrayID = as.numeric(args[1])

library(data.table)
library(Matrix)
library(xgboost)

#load environmental covariates:

EnvCovs = fread("~/ObesityCovOnlyCaucasiansIndepQC.txt",header = TRUE, colClasses = c("IID" = "character"))
#Select individuals from test data:

testIIDs = fread("~/TestDataIIDs.txt",header = TRUE, select = c("IID"), colClasses = c("IID" = "character"), sep = " ")

EnvCovs = subset(EnvCovs, EnvCovs$IID %in% testIIDs$IID)

EnvCovNames = colnames(EnvCovs)[2:ncol(EnvCovs)]

load("~/XGBmodels/XGB_FilteredIndepSNPsWithCovse0.05FromCVThres0.00084e0.05l1cl0.8ss0.8ct0.8m3Testdata4.RData")

model = RUN$Model

features = model$feature_names

snpnames = features[1:(length(features)-length(EnvCovNames))]

TestData = fread("~/RawFiles/TestDatae0.05Thres0.00085.raw", header = TRUE, sep = " ", drop = c("FID","PAT","MAT","SEX"), colClasses = c("IID" = "character"))

namesFromTestData = colnames(TestData)

#Find position in testdata where each snp from the model is positioned:

PositionOfsnpsInTestData = sapply(snpnames, function(x){ pos = grepl(x,namesFromTestData); ifelse(sum(pos) == 1, which(pos),0)})

snpsNotFound = snpnames[which(PositionOfsnpsInTestData == 0)]
#snpns not found is due to different allele frequencies and what is considered the minor allele.
#Try to find positions with the same snp-name, but different suffixes:
if(length(snpsNotFound)>0){
  snpsNotFoundClean = paste(sapply(strsplit(snpsNotFound,"_"),function(x) x[[1]]),"_",sep = "")
  samesnpname = sapply(snpsNotFoundClean, function(x) which(grepl(x,namesFromTestData)))
  PositionOfsnpsInTestData[which(PositionOfsnpsInTestData == 0)] = samesnpname
  
  colnames(TestData)[samesnpname] = snpsNotFound
}
PositionReordered = c(1,2,PositionOfsnpsInTestData)
#Reorder testdata such that the order of snpnames is the same as the order of snp names in TestData file:
TestData = TestData[,..PositionReordered]

#Merge environmental covariates to testdata:

TestData = merge(TestData,EnvCovs,by = "IID")
rm(EnvCovs)
rm(testIIDs)
gc()

TestData[,"IID":=NULL]

#Trace out the label (which is PHENOTYPE) (NB == 2 for HUNT data):
pheno_test = as.numeric(TestData$PHENOTYPE == 2)
TestData[,"PHENOTYPE":=NULL]

#Pick 100 people
inds = (arrayID*100-99):(arrayID*100)
pheno_test = pheno_test[inds]
TestData = TestData[inds,]

gc()
#XGBoost alogrithm requires the data structure to be a matrix:
TestDataAsMatrix = data.matrix(TestData)
rm(TestData)
gc()
#Construct the DMatrix:
