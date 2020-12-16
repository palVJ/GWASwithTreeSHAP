#Get all four matrices from all four models:
mats = list()
AllSNPs = c()
for(i in 1:4){
  
  #Sum all matrices of contributions for each interaction:
  library(Matrix)
  
  #load model:
  load(paste("~/RelInteractWithCovsIndepSNPsFromBestCVModelsTestData",i,"Array1.RData",sep = ""))
  AllSNPs = c(AllSNPs,rownames(RelInteractContribTest))
  mats[[i]] = RelInteractContribTest
  k = 1
  for(j in 2:470){
    
    if(file.exists(paste("~/RelInteractWithCovsIndepSNPsFromBestCVModelsTestData",i,"Array",j,".RData",sep = ""))){
      
      load(paste("~/RelInteractWithCovsIndepSNPsFromBestCVModelsTestData",i,"Array",j,".RData",sep = ""))
      mats[[i]] = mats[[i]] + RelInteractContribTest
      k = k + 1
    }
    
  }
  print(k)
}

AllSNPs = unique(AllSNPs)

#Create sparse matrix:

AverageRelInteractContrib = Matrix(data = 0, nrow = length(AllSNPs), ncol = length(AllSNPs), sparse = TRUE, dimnames = list(AllSNPs,AllSNPs))
HowManyTimesEvaluated = Matrix(data = 0, nrow = length(AllSNPs), ncol = length(AllSNPs), sparse = TRUE, dimnames = list(AllSNPs,AllSNPs))
for(i in 1:4){
  
  #find the pairs of non-zero interactions from matrix i:
  NonZeroElements = which(mats[[i]] != 0, arr.ind = TRUE)
  NonZeroElementsVectorIndices = which(mats[[i]] != 0)
  RowNames = rownames(mats[[i]])[which(mats[[i]] != 0, arr.ind = TRUE)[,1]]
  ColNames = colnames(mats[[i]])[which(mats[[i]] != 0, arr.ind = TRUE)[,2]]
  
  #Find the position of these interactions in the large matrix:
  PositionsInRow = sapply(RowNames, function(x) match(x,AllSNPs))
  PositionsInCol = sapply(ColNames, function(x) match(x,AllSNPs))
  
  #From array index to vector index in matrix:
  VectorIndex = sapply(1:length(PositionsInRow), function(x) (PositionsInCol[x]-1)*nrow(AverageRelInteractContrib)+PositionsInRow[x])
  
  #Insert interactions scores from matrix i to the large matrix:
  AverageRelInteractContrib[VectorIndex] = AverageRelInteractContrib[VectorIndex] + mats[[i]][NonZeroElementsVectorIndices]
  HowManyTimesEvaluated[VectorIndex] = HowManyTimesEvaluated[VectorIndex] + 1
}

#Multiply by the number of individuals that are investigated for each model:
HowManyTimesEvaluated = HowManyTimesEvaluated*(47000)

#For those interactions never evaluated add 1, as we want to divide AverageRelInteractContrib by HowManyTimesEvaluated and avoid dividing by zero.
HowManyTimesEvaluated[which(HowManyTimesEvaluated == 0)] = 1

AverageRelInteractContrib = AverageRelInteractContrib/HowManyTimesEvaluated

save(AverageRelInteractContrib,file = "~/AverageInteractEffects.RData")
