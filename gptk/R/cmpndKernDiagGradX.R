cmpndKernDiagGradX <-
function (kern, X) {
  X <- as.matrix(X)
  i <- 1
  funcName <- paste(kern$comp[[i]]$type, "KernDiagGradX", sep="")
  func <- get(funcName, mode="function")

  if ( !is.na(kern$comp[[i]]$index) ) {
    gX <- array(0, dim=dim(X))
    gX[,kern$comp[[i]]$index,] <- func(kern$comp[[i]], X[,kern$comp[[i]]$index])
  } else {
    gX <- func(kern$comp[[i]], X)
  }

  for ( i in seq(2, length=(length(kern$comp)-1)) ) {
    if ( !is.na(kern$comp[[i]]$index) ) {
      gX[,kern$comp[[i]]$index] <- gX[,kern$comp[[i]]$index] +  func(kern$comp[[i]], X[,kern$comp[[i]]$index])
    } else {
      gX <- gX + func(kern$comp[[i]], X)
    }
  }
 
  return (gX)
}
