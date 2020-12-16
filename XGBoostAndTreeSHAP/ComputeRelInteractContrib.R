#Sum all matrices of contributions for each interaction:
library(Matrix)

load("~/RelInteractWithCovsIndepSNPsFromBestCVModelsTestData3Array1.RData")
mat = RelInteractContribTest
k = 1
for(j in 2:470){
  
  if(file.exists(paste("~/RelInteractWithCovsIndepSNPsFromBestCVModelsTestData3Array",j,".RData",sep = ""))){
    
    load(paste("~/RelInteractWithCovsIndepSNPsFromBestCVModelsTestData3Array",j,".RData",sep = ""))
    mat = mat + RelInteractContribTest
    k = k + 1
  }
  
}

mat = mat/(k*100)
