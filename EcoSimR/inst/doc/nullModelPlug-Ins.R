## ----,eval=FALSE---------------------------------------------------------
#  foo <- function(speciesData,...){
#    ## Code to randomize here!
#  }

## ----,eval=FALSE---------------------------------------------------------
#  myAlgo <- function(speciesData) {
#  return(matrix(sample(speciesData), ncol=ncol(speciesData)))
#  }

## ----, eval = FALSE------------------------------------------------------
#  myAlgo2 <- function(speciesData,rowWeights,colWeights) {
#  matrixWeights <- outer(rowWeights,colWeights)
#  return(matrix(vector_sample(speciesData, weights =matrixWeights),ncol=ncol(speciesData)))
#  }

## ----,eval=FALSE---------------------------------------------------------
#  bar <- function(m){
#    return(ncol(unique(m,MARGIN=2)))
#  }

## ----, eval=FALSE--------------------------------------------------------
#  bar <- function(m, trim = FALSE){
#    if(trim){
#     m <- m[which(rowSums(m)>0),]
#     }
#    return(ncol(unique(m,MARGIN=2)))
#  }

## ------------------------------------------------------------------------
myAlgo <- function(speciesData) {
matrixWeights <- outer(rowSums(speciesData),rep(1,ncol(speciesData)))
return(matrix(vector_sample(speciesData, weights=matrixWeights),ncol=ncol(speciesData)))
}

## ----, message=FALSE,warning=FALSE,echo=FALSE----------------------------
library(EcoSimR)
## Simulate data
coocSimData <- ranMatGen(aBetaCol=0.5,bBetaCol=0.5,
aBetaRow=0.5,bBetaRow=0.5, numRows=30,numCols=30, mFill=0.25,abun=0,emptyRow=FALSE, emptyCol=FALSE)$m

coocOut <- null_model_engine(speciesData = coocSimData, algo = "myAlgo", metric = "checker",type="cooc",suppressProg = TRUE)


## ----,eval=FALSE---------------------------------------------------------
#  ## Simulate data
#  coocSimData <- ranMatGen(aBetaCol=0.5,bBetaCol=0.5, aBetaRow=0.5,bBetaRow=0.5, numRows=30,numCols=30,
#                            mFill=0.25,abun=0,emptyRow=FALSE, emptyCol=FALSE)$m
#  
#  coocOut <- null_model_engine(speciesData = coocSimData,nReps=1000,
#                               algo = "myAlgo", metric = "checker",type="cooc")
#  

## ----, fig.height=4,fig.width=4,fig.align='center'-----------------------
## Histogram plots
plot(coocOut,type='hist')

## ----,fig.height = 4,fig.width=6,fig.align='center'----------------------
## Matrix examples
plot(coocOut, type='cooc')

## Summary function
summary(coocOut)


## ----My algo and metric--------------------------------------------------

myAlgo <- function(speciesData,rowWeights,colWeights) {
matrixWeights <- outer(rowWeights,colWeights)
return(matrix(vector_sample(speciesData, weights =matrixWeights),ncol=ncol(speciesData)))
} 

## The C-Score
myMetric <- function(m) 
{
  m <- m[which(rowSums(m)>0),] # make calculation on submatrix with no missing species
  shared = tcrossprod(m)
  sums = rowSums(m)
  
  upper = upper.tri(shared)
  
  scores = (sums[row(shared)[upper]] - shared[upper])*
      (sums[col(shared)[upper]] - shared[upper])
  
  return(mean(scores))
}


## ----,message=FALSE,warning=FALSE,echo=FALSE-----------------------------

novelNull <- null_model_engine(coocSimData,algo="myAlgo",metric="myMetric",
                               algoOpts = list(rowWeights = runif(dim(coocSimData)[1]), 
                                               colWeights = runif(dim(coocSimData)[2])),suppressProg=T)


## ----,eval=FALSE---------------------------------------------------------
#  
#  novelNull <- null_model_engine(coocSimData,algo="myAlgo",metric="myMetric",
#                                 algoOpts = list(rowWeights = runif(dim(coocSimData)[1]),
#                                                 colWeights = runif(dim(coocSimData)[2])))

## ----, fig.height=4,fig.width=4,fig.align='center'-----------------------
summary(novelNull)
plot(novelNull)


