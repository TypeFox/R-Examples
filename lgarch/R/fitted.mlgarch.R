fitted.mlgarch <-
function(object, varma=FALSE, verbose=FALSE, ...)
{
  aux <- object$aux
  aux$verboseRecursion <- TRUE
  pars <- as.numeric(object$par.varma)
  mUhat <- mlgarchRecursion1(pars, aux)
  if(varma){
    result <- mUhat[,c(aux$m+1):c(2*aux$m)]
  }else{
    logsigma2 <- matrix(NA,aux$n,aux$m)
    sigma <- matrix(NA,aux$n,aux$m)
    i2 <- length(object$par)
    i1 <- i2 - aux$m + 1
    Elnz2 <- as.numeric(object$par[i1:i2])
    for(i in 1:aux$m){
      uhatnotzero <- as.numeric(mUhat[,i]!=0)
      lnz2adj <- mUhat[,i] + uhatnotzero*Elnz2[i]
      logsigma2[,i] <- mUhat[,c(i+aux$m)] - lnz2adj
      sigma[,i] <- exp(logsigma2[,i]/2)
    }
#    colnames(sigma) <- paste("sigmaFit",seq(1:aux$m),sep="")
    colnames(sigma) <- paste("sd",seq(1:aux$m),sep="")
    result <- sigma
    if(verbose){
#      colnames(logsigma2) <- paste("logsigma2Fit",seq(1,aux$m),sep="")
      colnames(logsigma2) <- paste("lnsd2no",seq(1,aux$m),sep="")
      result <- cbind(result,logsigma2)
    } #end if(verbose)
  } #end if(varma)else(..)
  result <- zoo(result, order.by=aux$y.index)
  return(result)
}
