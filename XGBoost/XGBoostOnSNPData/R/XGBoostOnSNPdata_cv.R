#' Use XGBoost to analyze SNP-data in a control/case GWAS setting using cross-validation in a HPC-environment
#' This function tries to find candidates of effects or interactions associated with a disease 
#' @param Trainingdata is a vector of .raw-SNP-data files, with indivdual per row, and one column showing
#' the phenotype of interest and the rest of the columns are the genotype for spesific SNPs along the genome.
#' @param covariatesTable is a data frame or data table of non-genetic covariates for each IID. 
#' @param CovariatesToInclude is a named vector of covariates (not including IID and SEX) to include.
#' @param Kfold is the number of folds wanted in cross-validation
#' @param features is a file consisting of one column of features names to be used from TrainData and ValData
#' @param hyperparameters is a list containing hyperparameters. See xgboost() package in R.
#' @param nrounds is the number of trees to be built during xgboost algorithm
#' @param early_stopping_rounds is the number of consecutive trees with no improvment in training error
#' before the algorithm stops to produce more trees. 
#' @return The mean test error after cross validation of xgboost-algorithm on SNP-data.
#' @export
#' @import data.table
#' @examples
#' XGBoostOnSNPdata_cv("~/directory/snpfile.txt",Kfold = 5, hyperparameters = list(eta = 0.3),
#' nrounds = 100, early_stopping_rounds = 10)



XGBoostOnSNPdata_cv = function(TrainData,ValData,features = NULL,covariatesTable = NULL,Kround,hyperparameters,nrounds, early_stopping_rounds,nthreads){
#XGBoost algorithm on one chromosome:
library(pryr)

#Get file for testdata:
TestDataFile = TrainData[Kround]
TrainingDataFiles = TrainData[-Kround]
#read Training data file which is a vector of files of subsets of the whole data to use for Training.
library(data.table)
if(is.null(features)){
subsetdata = lapply(TrainingDataFiles,function(x) fread(x, header = TRUE, sep = " ", drop = c("FID","PAT","MAT","SEX"),colClasses = c("IID" = "character")))
}
else if(features == "0"){
subsetdata = lapply(TrainingDataFiles,function(x) fread(x, header = TRUE, sep = " ", select = c("IID","PHENOTYPE"),colClasses = c("IID" = "character")))
}
else{

feat = as.character(read.table(file = features)[,1])
print(is.character(feat))
#SNPs may be in clean format (rs followed by number) or with suffix in addition to denote the minor allele:
if(grepl("_",feat[1])){

subsetdata = lapply(TrainingDataFiles,function(x) fread(x, header = TRUE, sep = " ", select = c("IID","PHENOTYPE",feat),colClasses = c("IID" = "character")))

}

else{

#read the names of the SNPs:
data = fread(TrainingDataFiles[1], header = TRUE, sep = " ", drop = c("FID","PAT","MAT","SEX","IID","PHENOTYPE"))
col_names = as.character(colnames(data))
rm(data)
gc()

#Find the position of the SNPs with same prefix as the ones in feat vector:
featWithSuffix = col_names[unlist(sapply(feat, function(x) which(grepl(paste(x,"_",sep = ""),col_names))))]

subsetdata = lapply(TrainingDataFiles,function(x) fread(x, header = TRUE, sep = " ", select = c("IID","PHENOTYPE",featWithSuffix),colClasses = c("IID" = "character")))

}

}
#merge data using rbindlist:
print(mem_used())
TrainingData = rbindlist(subsetdata)
rm(subsetdata)
gc()
print("The mem used after deleting Training data:")
#print(object_size(data))
print(mem_used())

#load covariatesfile if available:
if(is.data.frame(covariatesTable)){


#Find the covariates for training set:
TrainIIDs = TrainingData$IID
TrainCov = subset(covariatesTable,covariatesTable$IID %in% TrainIIDs)

#Remove participants due to kinship
#cov = subset(covariatesfile, !(as.character(unlist(cov[,"IID"])) %in% as.character(unlist(ParToRemove[,1]))))


#Merge SNP data and covariates data:
TrainingData = merge(TrainingData,TrainCov, by = "IID")
rm(TrainCov)
rm(TrainIIDs)
gc()
}
#The column IID is no longer needed
TrainingData[,c("IID"):=NULL]


#Create custom objective.

objective = function(logitpreds,dtrain){
labels = getinfo(dtrain, "label")
preds = 1/(1+exp(-logitpreds))

gradient = preds-labels
hessian = preds*(1-preds)

return(list(grad = gradient, hess = hessian))
}

#Create customized evaluation error.
evalerror = function(preds, data) {
  labels = getinfo(data, "label")
  #find rate of positives classified as negatives:
  typeIIerror = sum(labels == 1 & (preds>0)==0)/(sum(labels==1))
  #find type I error
  typeIerror = sum(labels == 0 & (preds>0) == 1)/(sum(labels==0))
  probs = 1/(1+exp(-preds))
  #Scores for cases:
  ScoresCases = probs[which(labels == 1)]
  ScoresControl = probs[which(labels == 0)]
  library(PRROC)
  PRAUC = pr.curve(ScoresCases,ScoresControl)$auc.integral
  return(list(metric = "PR-AUC",value = PRAUC, typeI = typeIerror, typeII = typeIIerror))
}


#Trace out the label (which is PHENOTYPE):
#set == 2 for HUNT data
pheno_train = as.numeric(TrainingData$PHENOTYPE == 2)

TrainingData[,"PHENOTYPE":=NULL]
#XGBoost alogrithm requires the data structure to be a matrix:
print("mem used before transforming data to matrix")
print(mem_used())
TrainingDataAsMatrix = data.matrix(TrainingData)
rm(TrainingData)
gc()
print("and after")
print(mem_used())

#Load the xgboost package:
library(xgboost)
#Construct the DMatrix:
dtrain = xgb.DMatrix(data=TrainingDataAsMatrix,label = pheno_train)

rm(TrainingDataAsMatrix)
gc()
#Download the test data:
print("mem used after making dtrain")
print(mem_used())
#TestData = rbindlist(list(Data_cases[-TrainingDataFromCases,],Data_controls[-TrainingDataFromControls,]))

if(is.null(features)){
TestData = fread(TestDataFile, header = TRUE, sep = " ", drop = c("FID","PAT","MAT","SEX"), colClasses = c("IID" = "character"))
}
else if(features == "0"){
TestData = fread(TestDataFile,header = TRUE,sep = " ", select = c("IID","PHENOTYPE"),colClasses = c("IID" = "character"))
}
else{

feat = as.character(read.table(features)[,1])

if(grepl("_",feat[1])){
TestData = fread(TestDataFile, header = TRUE, sep = " ", select = c("IID","PHENOTYPE",feat), colClasses = c("IID" = "character"))
}

else{

TestData = fread(TestDataFile, header = TRUE, sep = " ", select = c("IID","PHENOTYPE",featWithSuffix),colClasses = c("IID" = "character"))

}



}

if(is.data.frame(covariatesTable)){
#Include non-genetic covariates:

TestIIDs = TestData$IID
TestCov = subset(covariatesTable, covariatesTable$IID %in% TestIIDs)

#merge data:
TestData = merge(TestData,TestCov, by = "IID")

rm(TestCov)
rm(TestIIDs)
gc()

}
#IID no longer needed:
TestData[,c("IID"):=NULL]

#Trace out the label (which is PHENOTYPE) (NB == 2 for HUNT data):
pheno_test = as.numeric(TestData$PHENOTYPE == 2)
TestData[,"PHENOTYPE":=NULL]
#XGBoost alogrithm requires the data structure to be a matrix:
TestDataAsMatrix = data.matrix(TestData)
rm(TestData)
gc()
#Construct the DMatrix:
dtest = xgb.DMatrix(data=TestDataAsMatrix,label = pheno_test)
rm(TestDataAsMatrix)
gc()
#Train the model:
gc()
print("mem used just before starting model")
print(mem_used())
###############TRAIN MODEL HERE######################
model =  xgb.train(data = dtrain, obj = objective,feval = evalerror,nrounds = nrounds,
params = hyperparameters,
nthread = nthreads,maximize = TRUE,watchlist = list(test=dtest),verbose = 3, early_stopping_rounds = early_stopping_rounds)
#models[[i]] = model
#########################################################


#FEATURE IMPORTANCE: Find which covariates are the most important in each training:
ImportanceRoundi = xgb.importance(model = model)
#ImportanceRoundi[,"Round"] = rep(Kround,nrow(ImportanceRoundi))
#Importance = rbind(Importance,ImportanceRoundi)


rm(dtrain)
gc()
# generate predictions for our held-out validation data:
if(is.null(features)){
ValidData = fread(ValData, header = TRUE, sep = " ", drop = c("FID","PAT","MAT","SEX"), colClasses = c("IID"="character"))
}
else if(features == "0"){
ValidData = fread(ValData,header = TRUE,sep = " ", select = c("IID","PHENOTYPE"),colClasses = c("IID" = "character"))
}
else{

feat = as.character(read.table(features)[,1])

if(grepl("_",feat[1])){
ValidData = fread(ValData, header = TRUE, sep = " ", select = c("IID","PHENOTYPE",feat), colClasses = c("IID"="character"))
}

else{

ValidData = fread(ValData, header = TRUE, sep = " ", select = c("IID","PHENOTYPE",featWithSuffix),colClasses = c("IID" = "character"))

}


}
if(is.data.frame(covariatesTable)){
#Include non-genetic covariates:

ValidIIDs = ValidData$IID
ValidCov = subset(covariatesTable, covariatesTable$IID %in% ValidIIDs)

#merge data:
ValidData = merge(ValidData,ValidCov, by = "IID")

rm(ValidCov)
rm(ValidIIDs)
gc()

}
#IID no longer needed:
ValidData[,c("IID"):=NULL]


#Trace out the label (which is PHENOTYPE) (NB == 2 for HUNT data):
pheno_val = as.numeric(ValidData$PHENOTYPE == 2)
ValidData[,"PHENOTYPE":=NULL]
gc()
#XGBoost alogrithm requires the data structure to be a matrix:
ValDataAsMatrix = data.matrix(ValidData)
rm(ValidData)
gc()
#Construct the DMatrix:
dVal = xgb.DMatrix(data=ValDataAsMatrix,label = pheno_val)
rm(ValDataAsMatrix)
gc()

pred_val = predict(model, dVal)



#compute validation error:
evaluate = evalerror(pred_val,dVal)
Valerror = evaluate$value
typeIerror = evaluate$typeI
typeIIerror = evaluate$typeII

rm(dVal)
print(mem_used())
gc()
#}
  

return(list(Score = Valerror,TypeIerror = typeIerror, TypeIIerror = typeIIerror,Importance = ImportanceRoundi, Model = model))
}
