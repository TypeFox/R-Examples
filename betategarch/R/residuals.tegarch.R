residuals.tegarch <-
function(object, standardised=TRUE, ...)
{
y.index <- index(object$y)
y <- coredata(object$y)

if(object$model["components"]==1){
  if(object$model["skew"]==0){
    pars <- c(object$par, 1)
  }else{ pars <- object$par }
  if(length(pars) < 6){
    pars <- c(pars[1:3], 0, pars[4:5])
  }
  out <- tegarchRecursion(y, omega=pars[1], phi1=pars[2],
    kappa1=pars[3], kappastar=pars[4], df=pars[5],
    skew=pars[6], lambda.initial=object$lambda.initial,
    c.code=TRUE, verbose=TRUE, aux=NULL)
  out <- zoo(out, order.by=y.index)
}else{
  if(object$model["skew"]==0){
    pars <- c(object$par, 1)
  }else{ pars <- object$par }
  out <- tegarchRecursion2(y, omega=pars[1], phi1=pars[2],
    phi2=pars[3], kappa1=pars[4], kappa2=pars[5],
    kappastar=pars[6], df=pars[7], skew=pars[8],
    lambda.initial=object$lambda.initial, c.code=TRUE,
    verbose=TRUE, aux=NULL)
  out <- zoo(out, order.by=y.index)
}
if(standardised){
  out <- out[,"residstd"]
}else{
  out <- out[,"epsilon"]
}
return(out)
}
