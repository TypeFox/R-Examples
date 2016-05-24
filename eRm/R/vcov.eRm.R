`vcov.eRm` <-
function(object,...) 
{
  if (any(is.na(object$se.eta))) {
    vcmat <- NA 
  } else {
    vcmat <- (solve(object$hessian))      #VC-matrix of the parameter estimates
  }
  return(vcmat)
}

