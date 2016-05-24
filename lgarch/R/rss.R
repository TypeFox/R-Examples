rss <-
function(object, ...)
{
  if(class(object)!="lgarch" && class(object)!="armax")
    stop("object does not belong to the lgarch class")
  
  if(class(object)=="lgarch"){
    pars <- as.numeric(object$par.arma)
    if(object$aux$mean.correction){ pars[1] <- 0 }
    uhat <- lgarchRecursion1(pars, object$aux)
    if(object$aux$yzeron > 0){
      uhat <- uhat[-object$aux$yzerowhere]
    }
  }
  
  return(sum(uhat^2))
}
