#  File R/logLik.stergm.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################

logLik.stergm<-function(object, add=FALSE, force.reeval=FALSE, eval.loglik=add || force.reeval, control=control.logLik.stergm(), ...){
  check.control.class()
  if(object$estimate=="EGMME") stop("Log-likelihood for ",object$estimate," is not meaningful.")
  
  if(!is.null(control$seed))  set.seed(as.integer(control$seed))

  if(add){
    object$formation.fit <- logLik(object$formation.fit, add=add, force.reeval=force.reeval, eval.loglik = add || force.reeval, control=control$control.form)
    object$dissolution.fit <- logLik(object$dissolution.fit, add=add, force.reeval=force.reeval, eval.loglik = add || force.reeval, control=control$control.diss)
    
    object
  }else{
    llk.form <- logLik(object$formation.fit, add=add, force.reeval=force.reeval, eval.loglik = eval.loglik, control=control$control.form)
    llk.diss <- logLik(object$dissolution.fit, add=add, force.reeval=force.reeval, eval.loglik = eval.loglik, control=control$control.diss)

    llk <- llk.form + llk.diss
    class(llk) <- "logLik"
    attr(llk,"df") <- attr(llk.form,"df") + attr(llk.diss,"df")
    attr(llk,"nobs") <- attr(llk.form,"nobs") + attr(llk.diss,"nobs")
    
    llk
  }
}

logLikNull.stergm <- function(object, control=control.logLik.stergm(), ...){
    llk.form <- logLikNull(object$formation.fit, control=control$control.form)
    llk.diss <- logLikNull(object$dissolution.fit, control=control$control.diss)

    llk <- llk.form + llk.diss
    class(llk) <- "logLik"
    attr(llk,"df") <- attr(llk.form,"df") + attr(llk.diss,"df")
    attr(llk,"nobs") <- attr(llk.form,"nobs") + attr(llk.diss,"nobs")
    
    llk
}
