cmpndKernGradX <-
function (kern, X, X2) {
  i <- 1
  funcName <- paste(kern$comp[[i]]$type, "KernGradX", sep="")
  func <- get(funcName, mode="function")

  if ( !is.na(kern$comp[[i]]$index) ) {
    gX <- array(0, dim=c(dim(as.array(X2))[1], dim(as.array(X2))[1],
                     dim(as.array(X))[1]))
    gX[,kern$comp[[i]]$index,] <- func(kern$comp[[i]], X[,kern$comp[[i]]$index], X2[,kern$comp[[i]]$index])
  } else {
    gX <- func(kern$comp[[i]], X, X2)
  }

  for ( i in seq(2, length=(length(kern$comp)-1)) ) {
    funcName <- paste(kern$comp[[i]]$type, "KernGradX", sep="")
    func <- get(funcName, mode="function")
    if ( !is.na(kern$comp[[i]]$index) ) {
      gX[,kern$comp[[i]]$index,] <- gX[,kern$comp[[i]]$index,] +  func(kern$comp[[i]], X[,kern$comp[[i]]$index], X2[,kern$comp[[i]]$index])
    } else {
      gX <- gX + func(kern$comp[[i]], X, X2)
    }
  }
 
  return (gX)
}
