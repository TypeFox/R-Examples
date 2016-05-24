fitted.lgarch <-
function(object, verbose=FALSE, ...)
{
  object$aux$verboseRecursion <- TRUE
  pars <- as.numeric(object$par.arma)
  if(object$aux$mean.correction){ pars[1] <- 0 }
  mUhat <- lgarchRecursion1(pars, object$aux)
  lnz2adj <- mUhat[,1] + object$par["Elnz2"]
  if(object$aux$mean.correction){
    logsigma2 <- mUhat[,2] + object$aux$Elny2 - lnz2adj
  }else{
    logsigma2 <- mUhat[,2] - lnz2adj
  }
  sigma <- exp(logsigma2/2)
  if(verbose){
    zhat <- object$aux$y/sigma
    result <- cbind(sigma, logsigma2, zhat, mUhat[,1])
    colnames(result) <- c("sd", "lnsd2", "zhat", "uhat")
  }else{
    result <- sigma
  }
  result <- zoo(result, order.by=object$aux$y.index)
  return(result)
}
