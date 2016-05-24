deviance.nelm <-
function(object, ...)
{
  nme  <- length(object$erg_distr$mean_est) - 1
  nva  <- object$ctrl$sigmaest *(length(object$erg_distr$sig_est) -1)
  npar <- ncol(object$reshOBJ$Qmat) + nme + nva - length(object$ctrl$Clist)
  
structure(2*object$last_mstep$value, df=npar)
}



## logLik ##

logLik.nelm <- 
  function(object,...)
  {
    
    ### number of parameters
    nme  <- length(object$erg_distr$mean_est) - 1
    nva  <- object$ctrl$sigmaest * (length(object$erg_distr$sig_est) -1)
    npar <- ncol(object$reshOBJ$Qmat) + nme + nva  - length(object$ctrl$Clist)
    # number of observations
    nobs <- sum(sapply(object$reshOBJ$d,nrow))
    
    return(structure((-1)*object$last_mstep$value, df=npar, nobs=nobs))
    
  }


