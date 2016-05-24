residuals.lgarch <-
function(object, arma=FALSE, ...)
{
  if(arma){
    pars <- as.numeric(object$par.arma)
    if(object$aux$mean.correction){ pars[1] <- 0 }
    result <- lgarchRecursion1(pars, object$aux)
  }else{
    sigma <- coredata(fitted.lgarch(object))
    result <- object$aux$y/sigma
  }
  result <- zoo(result, order.by=object$aux$y.index)
  return(result)
}
