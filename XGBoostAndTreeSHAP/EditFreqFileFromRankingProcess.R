#load freq-file:

freqfile = read.table("~/AllSNPsUsedAfterRankingWithCovsIndepSNPs.frq",header = TRUE)

#load Scores of features file:

scorefile = read.table("~/ScoreFilteredWithCovsIndepSNPs.txt",header = TRUE)


#Adjust the freq file such that the column of "MAF" in .frq-file is now a measure of importance based on the score in scorefile for each SNP.
#The larger score, the larger "MAF". Simply set MAF equal to score:


#For each snp in freqfile, change the MAF value to its corresponding score:

for(i in 1:nrow(freqfile)){

#Find the corresponding snp in the scorefile and swap the MAF-value with score:
snp = which(grepl(paste(as.character(freqfile$SNP[i]),"_",sep = ""),as.character(scorefile$Feature)))

#The snp may be named with two suffixes depending on which allele became minor allele in traing data (may occur for MAF close to 0.5).
#Set the score equal to the average when this happens:
if(length(snp)>1){
print(scorefile$Feature[snp])
}
freqfile$MAF[i] = sum(scorefile$Score[snp])/length(snp)


}


write.table(freqfile,"~/AllSNPsUsedAfterCVWithCovsIndepSNPAdjusted.frq",col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
