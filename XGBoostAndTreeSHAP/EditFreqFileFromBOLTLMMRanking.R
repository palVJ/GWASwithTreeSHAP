#load freq-file:

freqfile = read.table("~/FreqFileAllSNPs.frq",header = TRUE)

#load Scores of features file:
library(data.table)
scorefile = fread("~/BOLT_output_obesity_withcovariates",header = TRUE,
 colClasses = c("SNP" = "character","P_BOLT_LMM_INF" = "numeric"), select = c("SNP","P_BOLT_LMM_INF"))

#Only look at top 40 000 SNPs:
scorefile = scorefile[1:40000,]

#get the reduced freqfile with SNPs in the scorelist:
freqfile = subset(freqfile, freqfile$SNP %in% scorefile$SNP)


#Adjust the freq file such that the column of "MAF" in .frq-file is now a measure of importance based on the score in scorefile for each SNP.
#The larger score, the larger "MAF". Simply set MAF equal to score:

#First change scorefile such that scores are between 0 and 0.5 instead of p-values:
seq = a = (0:(nrow(scorefile)-1))/1000
scorefile$P_BOLT_LMM_INF = 0.5*exp(-seq)
#For each snp in freqfile, change the MAF value to its corresponding score:

for(i in 1:nrow(freqfile)){

#Find the corresponding snp in the scorefile and swap the MAF-value with score:
snp = which(as.character(freqfile$SNP[i]) == as.character(scorefile$SNP))

freqfile$MAF[i] = scorefile$P_BOLT_LMM_INF[snp]


}


write.table(freqfile,"~/AllSNPsUsedAfterBOLTLMMAdjusted.frq",col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
