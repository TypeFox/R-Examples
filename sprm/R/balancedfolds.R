balancedfolds <- function (y, nfolds) 
{
  nA <- sum(y==1)
  nB <- sum(y==-1)
  foldsA <- cvFolds(nA, K = nfolds, R = 1, type = "random")
  foldsB <- cvFolds(nB, K = nfolds, R = 1, type = "random")
  
  folds <- list()
  for(f in 1:nfolds){
    folds[[f]] <- c(which(y==1)[foldsA$subsets[foldsA$which==f,]], which(y==-1)[foldsB$subsets[foldsB$which==f,]])
  }
  return(folds)
}

