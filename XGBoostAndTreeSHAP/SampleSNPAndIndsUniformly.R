#This script is used to sample which snps and individuals to include in XGBoost model:
args <- commandArgs(trailingOnly = TRUE)
#number of snps to sample:

N = as.numeric(args[1])
shID = as.numeric(args[3])
#Read the name of all the snps:
library(data.table)

snps = fread("~/RawFiles/MAFfreqObesity.frq.cc",select = c("SNP"),
colClasses = c("SNP"="character"), header = TRUE)

#Sample N snps:

Subsample = snps$SNP[sample(x=nrow(snps),size=N,replace = FALSE)]

#write table:
table = data.frame(Subsample)
fwrite(table,file = paste("~/RawFiles/snpsToIncludeBeforeIndep",shID,".txt",sep = ""),row.names = FALSE,quote = FALSE,
col.names = FALSE)


#Read the phenotypefile:
ids =  fread("~/PhenoFileOnlyCaucasiansIndsQC.txt",colClasses = c("FID"="character","IID" = "character"),sep = " ",header = TRUE)

#Read RankingData IIDs:

RankingDataIIDs = fread("~/RankingDataIIDs.txt",header = TRUE,select = c("IID"), colClasses = c("IID" = "character"))$IID

ids = subset(ids, ids$IID %in% RankingDataIIDs)

#Number of individuals to include:
S = as.numeric(args[2])
#Proportion of cases:
IndsInTotal = sum(!is.na(ids$Obese))
CCR = sum(ids$Obese,na.rm = T)/IndsInTotal
#Number of cases out of the total S individuals to include in order to preserve the proportion of cases:
NrOfCases = round(S*CCR)

#Sample randomly controls from population:
Controls = which(ids$Obese == 0)
ControlsToInclude = sample(x = Controls, size = S-NrOfCases, replace = FALSE)
Cases = which(ids$Obese == 1)
CasesToInclude = sample(x = Cases,size = NrOfCases, replace = FALSE)

IIDsToInclude = ids$FID[c(ControlsToInclude,CasesToInclude)]

IIDsToInclude = data.frame(FID = IIDsToInclude, IID = IIDsToInclude)
fwrite(IIDsToInclude,file = paste("~/RawFiles/IIDsToIncludeIndepSNPs",shID,".txt",sep = ""),row.names = FALSE,quote = FALSE,
col.names = TRUE, sep = "\t")
