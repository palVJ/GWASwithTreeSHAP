#Split data in stratified subsets for use to create training data, validation data and test data:
args <- commandArgs(trailingOnly = TRUE)
#load data:
shID = args[2]
library(data.table)


data = fread(paste("~/RawFiles/SNP_additive_data_unifIndepSNPs",shID,".raw",sep = ""),header = TRUE, sep = " ",nThread=3)

#The number of folds to create:
folds = as.numeric(args[1])

Cases = which(data$PHENOTYPE == 2)
Controls = which(data$PHENOTYPE == 1)
NrOfCases = length(Cases)
NrOfControls = length(Controls)
#Divide in K folds for cross validation:

NrOfSamplesPerFold = floor(nrow(data)/folds)
#Amount of cases in each fold such that the true proportion of cases is preserved in each fold:
CasesPerFold = floor(NrOfCases/folds)
ControlsPerFold = floor(NrOfControls/folds)

#Randomly sample from cases and controls:
RandSampleFromCases = sample(Cases,NrOfCases,replace = FALSE)
RandSampleFromControls =  sample(Controls,NrOfControls,replace = FALSE)
#Load the samples (rownumbers) for each fold in a list:
SamplesPerFold = list()
temp1 = 1
temp2 = 1
for (i in 1:(folds-1)){
SamplesPerFold[[i]] = list(RandSampleFromCases[temp1:(temp1+CasesPerFold-1)],RandSampleFromControls[temp2:(temp2+ControlsPerFold-1)])
temp1 = temp1 + CasesPerFold
temp2 = temp2 + ControlsPerFold

#Slice the data with rows equal to those in SamplesPerFold[[i]]
subsetdata = data[unlist(SamplesPerFold[[i]])]

#Write table:
fwrite(subsetdata,file = paste("~/RawFiles/","SNP_add_SubsetIndepSNPs",shID,"_Trai",i,".raw",sep = ""),quote=FALSE,row.names=FALSE,sep=" ",na="NA",nThread=3)
rm(subsetdata)
gc()
}
gc()
#The validation data:
SamplesPerFold[[folds]] = list(RandSampleFromCases[temp1:length(RandSampleFromCases)],RandSampleFromControls[temp2:length(RandSampleFromControls)])
subsetdata = data[unlist(SamplesPerFold[[folds]])]
fwrite(subsetdata,file = paste("~/RawFiles/","SNP_add_SubsetIndepSNPs",shID,"_ValData",".raw",sep = ""),quote=FALSE,row.names=FALSE,sep=" ",na="NA",nThread=3)
