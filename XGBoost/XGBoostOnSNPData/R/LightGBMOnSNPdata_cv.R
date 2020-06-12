#' Use LightGBM to analyze SNP-data in a control/case GWAS setting using cross-validation in a HPC-environment
#' This function tries to find candidates of effects or interactions associated with a disease
#' @param Trainingdata is a vector of .raw-SNP-data files, with indivdual per row, and one column showing
#' the phenotype of interest and the rest of the columns are the genotype for spesific SNPs along the genome.
#' @param covariatesfile is a table of non-genetic (PCs are OK) covariates for each IID.
#' @param CovariatesToInclude is a named vector of covariates (not including IID and SEX) to include.
#' @param Kfold is the number of folds wanted in cross-validation
#' @param hyperparameters is a list containing hyperparameters. See xgboost() package in R.
#' @param nrounds is the number of trees to be built during xgboost algorithm
#' @param early_stopping_rounds is the number of consecutive trees with no improvment in training error
#' before the algorithm stops to produce more trees.
#' @return The mean test error after cross validation of xgboost-algorithm on SNP-data.
#' @export
#' @import data.table
#' @examples
#' LightGBMOnSNPdata_cv("~/directory/snpfile.txt",Kfold = 5, hyperparameters = list(eta = 0.3),
#' nrounds = 100, early_stopping_rounds = 10)



LightGBMOnSNPdata_cv = function(TrainData,ValData,covariatesfile = NULL,CovariatesToInclude = NULL,Kround,hyperparameters,nrounds, early_stopping_rounds,nthreads,max_depth){
#XGBoost algorithm on one chromosome:
library(pryr)

#Get file for testdata:
TestDataFile = TrainData[Kround]
TrainingDataFiles = TrainData[-Kround]
#read Training data file which is a vector of files of subsets of the whole data to use for Training.
library(data.table)
subsetdata = lapply(TrainingDataFiles,function(x) fread(x, header = TRUE, sep = " ", drop = c("IID","FID","PAT","MAT","SEX")))
#merge data using rbindlist:
print(mem_used())
TrainingData = rbindlist(subsetdata)
rm(subsetdata)
gc()
print("The mem used after deleting Training data:")
#print(object_size(data))
print(mem_used())
#Remove individuals where phenotype is not known
#data = subset(data,as.vector(unlist(data[,"PHENOTYPE"]))!=-9)
#IIDsLeft = as.vector(unlist(data[,"IID"]))

#download IDs of participants to remove due to kinship:

#ParToRemove = fread("~/CleaningTheData/ParticipantsToRemoveDueToKinship.txt",sep = "\t",colClasses = c("character","NULL"))
"$shID"
#Remove participants due to kinship from SNP data
#data = subset(data, !(as.character(unlist(data[,"IID"])) %in% ParToRemove[,1]))

#Remove columns that are not of any interest:
#Trainingdata[,c("FID","PAT","MAT","SEX"):=NULL]

#load covariatesfile if available:
if(!is.null(nrow(covariatesfile))){
cov = fread(covariatesfile,header = TRUE,colClasses = c(IID = "character"), sep = "\t", select = c("IID",CovariatesToInclude))

#Only include IIDs with known phenotype
cov = subset(cov,cov$IID %in% IIDsLeft)

#Remove participants due to kinship
#cov = subset(covariatesfile, !(as.character(unlist(cov[,"IID"])) %in% as.character(unlist(ParToRemove[,1]))))


#Merge SNP data and covariates data:

Trainingdata = merge(data,cov, by = "IID")
#The column IID is no longer needed and can not be used as covariate:
}
#Trainingdata[,c("IID"):=NULL]

#Split into cases and control:
#SET == 2 FOR HUNT DATA in SNP data file
#Cases = which(data$PHENOTYPE == 2)
#Data_cases = data[Cases,]
#Data_controls = data[-Cases,]
#NrOfCases = nrow(Data_cases)
#NrOfControls = nrow(Data_controls)
#Divide in k folds for cross validation:
#folds = Kfold
#NrOfSamplesPerFold = floor(nrow(data)/folds)
#Amount of cases in each fold such that the true proportion of cases is preserved in each fold:
#CasesPerFold = floor(NrOfCases/folds)
#ControlsPerFold = floor(NrOfControls/folds)
#rm(Trainingdata)
#gc()
#Randomly sample from cases and controls:
#RandSampleFromCases = sample(1:NrOfCases,NrOfCases,replace = FALSE)
#RandSampleFromControls =  sample(1:NrOfControls,NrOfControls,replace = FALSE)
#Load the samples (rownumbers) for each fold in a list:
#SamplesPerFold = list()
#temp1 = 1
#temp2 = 1
#for (i in 1:(folds-1)){
#SamplesPerFold[[i]] = list(RandSampleFromCases[temp1:(temp1+CasesPerFold-1)],RandSampleFromControls[temp2:(temp2+ControlsPerFold-1)])
#temp1 = temp1 + CasesPerFold
#temp2 = temp2 + ControlsPerFold
#}

#SamplesPerFold[[folds]] = list(RandSampleFromCases[temp1:length(RandSampleFromCases)],RandSampleFromControls[temp2:length(RandSampleFromControls)])


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
  return(list(name = "PR-AUC",value = PRAUC,typeI = typeIerror,typeII = typeIIerror, higher_better = TRUE))
}


#For each loop, compute test error. The mean test error will be a measure
#of how well the model performs:
#testerror = c()
#typeIerrors = c()
#typeIIerrors = c()
#Importance = data.frame(matrix(ncol = 5, nrow = 0))
#colnames(Importance) = c("Round","Feature","Gain","Cover","Frequency")

#Save each trained lightgbm model in the list models

#models = list()
#for (i in 1:folds){
#Fold i is decided to be test sample, in other words merge data from the rest of the folds:
#TrainingDataFromControls = (1:NrOfControls)[-unlist(SamplesPerFold[[i]][2])]
#TrainingDataFromCases = (1:NrOfCases)[-unlist(SamplesPerFold[[i]][1])]
#Bind training data from cases and control together:
#TrainingData = rbindlist(list(Data_cases[TrainingDataFromCases,],Data_controls[TrainingDataFromControls,]))
#Trace out the label (which is PHENOTYPE):
#set == 2 for HUNT data
pheno_train = as.numeric(TrainingData$PHENOTYPE == 2)

TrainingData[,"PHENOTYPE":=NULL]
#XGBoost alogrithm requires the data structure to be a matrix:
print("mem used before transforming data to matrix")
print(mem_used())
#TrainingDataAsMatrix = data.matrix(TrainingData)
#rm(TrainingData)
gc()
print("and after")
print(mem_used())

#Load the lightgbm package:
library(lightgbm)
TrainingDataAsMatrix= as.matrix(lgb.prepare2(TrainingData))
rm(TrainingData)
gc()
#Construct the lgb.Dataset:
dtrain = lgb.Dataset(data=TrainingDataAsMatrix,label = pheno_train)

rm(TrainingDataAsMatrix)
gc()
#Download the test data:
print("mem used after making dtrain")
print(mem_used())
#TestData = rbindlist(list(Data_cases[-TrainingDataFromCases,],Data_controls[-TrainingDataFromControls,]))
TestData = fread(TestDataFile, header = TRUE, sep = " ", drop = c("FID","IID","PAT","MAT","SEX"))

#Trace out the label (which is PHENOTYPE) (NB == 2 for HUNT data):
pheno_test = as.numeric(TestData$PHENOTYPE == 2)
TestData[,"PHENOTYPE":=NULL]
#LightGBM alogrithm requires the data structure to be a matrix:
#TestDataAsMatrix = data.matrix(TestData)
#rm(TestData)
gc()
#Construct the lgb.Dataset:
TestDataAsMatrix= as.matrix(lgb.prepare2(TestData))
rm(TestData)
gc()
dtest = lgb.Dataset(data=TestDataAsMatrix,label = pheno_test)
rm(TestDataAsMatrix)
gc()
#Train the model:
gc()
print("mem used just before starting model")
print(mem_used())
###############TRAIN MODEL HERE######################
model =  lgb.train(data = dtrain, obj = objective,eval = evalerror,nrounds = nrounds,
params = hyperparameters,
num_threads = nthreads,valids = list(test=dtest),verbose = 1, early_stopping_rounds = early_stopping_rounds,max_depth = max_depth)
#models[[i]] = model
#########################################################


#FEATURE IMPORTANCE: Find which covariates are the most important in each training:
ImportanceRoundi = lgb.importance(model = model)
#ImportanceRoundi[,"Round"] = rep(Kround,nrow(ImportanceRoundi))
#Importance = rbind(Importance,ImportanceRoundi)


rm(dtrain)
gc()
# generate predictions for our held-out validation data:
ValidData = fread(ValData, header = TRUE, sep = " ", drop = c("FID","IID","PAT","MAT","SEX"))

#Trace out the label (which is PHENOTYPE) (NB == 2 for HUNT data):
pheno_val = as.numeric(ValidData$PHENOTYPE == 2)
ValidData[,"PHENOTYPE":=NULL]
gc()
#LightGBM alogrithm requires the data structure to be a matrix:
ValDataAsMatrix = data.matrix(ValidData)
rm(ValidData)
gc()
#Construct the lgb.Dataset:
dVal = lgb.Dataset(data=ValDataAsMatrix,label = pheno_val)
#rm(ValDataAsMatrix)
#gc()

pred_val = predict(model, ValDataAsMatrix,rawscore = TRUE)
rm(ValDataAsMatrix)
gc()
#compute validation error:
evaluate = evalerror(pred_val,dVal)
Valerror = evaluate$value
typeIerror = evaluate$typeI
typeIIerror = evaluate$typeII

print(mem_used())
gc()
#}


return(list(Score = Valerror,TypeIerror = typeIerror,TypeIIerror = typeIIerror,Importance = ImportanceRoundi, Model = model))
}

