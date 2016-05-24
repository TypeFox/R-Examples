

 #######################
 #### summary.ldbglm ####
 #######################

 ## Description:
 ##     Summary of an object of class ldbglm.
 ##         - calculate the pearson residuals
 ##         - the deviance residuals
 ##         - the dispersion
 ##         - the hat trace
 ##        

summary.ldbglm <-function(object,dispersion=NULL,...){

    if (!inherits(object, "ldbglm")) 
      stop("use only with \"ldbglm\" objects")
 
    # recover attributes rdf, weights, residuals, etc. of dbglm object
    z <- object
    nobs <- length(z$y)
    weights <- object$weights

    y <-object$y
    mu <-object$fitted.values
    G<-attr(object,"G")

    # trace_hat
    t_hat <- sum(diag(z$S))

    # residual sigma
    sigma<-attr(z,"sigma")

    # matched call
    call <-z$call

    # family
    family <-object$family
    dev.resids <-family$dev.resids
    mu.eta <-family$mu.eta
   
    kind.kernel <- switch(attr(object,"kind.of.kernel"),
     "(1) Epanechnikov",
	   "(2) Biweight",
		 "(3) Triweight",
		 "(4) Normal",
		 "(5) Triangular",
		 "(6) Uniform")
    
      crit.value <- NULL
      if (!is.null(attr(object,"OCV_opt")))
        crit.value <- attr(object,"OCV_opt") 
      if (!is.null(attr(object,"GCV_opt")))
        crit.value <- attr(object,"GCV_opt") 
      if (!is.null(attr(object,"AIC_opt")))
        crit.value <- attr(object,"AIC_opt") 
      if (!is.null(attr(object,"BIC_opt")))
        crit.value <- attr(object,"BIC_opt")
         
   ## definitions analogous to those in summary.dbglm
   # residuals, residuals degree of freedom, null deviance 
   wtdmu         <-sum(weights*y)/sum(weights)
   n.ok          <-nobs-sum(weights==0)
   null.deviance <-sum(dev.resids(y,wtdmu,weights))
   df.null       <-n.ok-1
   df.residual   <-n.ok-t_hat
   residual.deviance <- sum(dev.resids(y,mu,weights))


   # pearson residuals
   pears.resid<-(y-mu)/sqrt(family$variance(mu))
   
   # deviance residuals 
   deviance.resid<-sign(y-mu)*sqrt(dev.resids(y,mu,weights))

   # residuals for defining "dispersion"
   eta<-family$linkfun(mu)
   residuals <- (y-mu)/mu.eta(eta)
      
   # dispersion
    est.disp <- FALSE
    if (is.null(dispersion)){
     dispersion <- if (family$family %in% c("poisson","binomial")) 1
                  else if (df.residual > 0) {
                    est.disp <- TRUE
                    if (any(weights == 0))
                      warning("observations with zero weight not used for calculating dispersion")
                    sum((weights * residuals^2)[weights > 0])/df.residual
                 }
                 else {
                   est.disp <- TRUE
                   NaN
                 }
   } else{
    if (!is.numeric)
      stop("dispersion must be numeric")
   }
   ## end of definitions analogous to those in summary.dbglm   
    summary.dist1 <- summary(object$dist1) 
    
    # list to be returned
    ans <- list(nobs=nobs,trace.hat=t_hat,call=call, family=family,y=z$y,
                residual.deviance=residual.deviance, df.residual=df.residual,
                null.deviance=null.deviance, df.null=df.null,
                deviance.resid=deviance.resid, pears.resid=pears.resid,
                dispersion=dispersion, residuals=residuals,
                kind.kernel=kind.kernel,method.h=attr(object,"method.h"),
                h.opt=object$h.opt, crit.value=crit.value,
                summary.dist1=summary(as.numeric(object$dist1^.5)),
                percentile.h.opt=100*mean(object$dist1^.5<=object$h.opt))
 
   
    class(ans)<-"summary.ldbglm"

    return(ans)
}



