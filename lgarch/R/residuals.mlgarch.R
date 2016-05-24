residuals.mlgarch <-
function(object, varma=FALSE, ...)
{
  aux <- object$aux
  aux$verboseRecursion <- TRUE
  pars <- as.numeric(object$par.varma)
  mUhat <- mlgarchRecursion1(pars, aux)
  if(varma){
    result <- mUhat[,1:aux$m]
  }else{
    result <- matrix(NA,aux$n,aux$m)
    i2 <- length(object$par)
    i1 <- i2 - aux$m + 1
    Elnz2 <- as.numeric(object$par[i1:i2])
    for(i in 1:aux$m){
      uhatnotzero <- as.numeric(mUhat[,i]!=0)
      lnz2adj <- mUhat[,i] + uhatnotzero*Elnz2[i]
      logsigma2 <- mUhat[,c(i+aux$m)] - lnz2adj
      sigma <- exp(logsigma2/2)
      result[,i] <- aux$y[,i]/sigma
    }
    colnames(result) <- paste("z",seq(1,aux$m),sep="")
  }
  result <- zoo(result, order.by=aux$y.index)
  return(result)
}
