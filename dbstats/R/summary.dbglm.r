
 ####################### 
 #### summary.dbglm ####
 #######################

 ## Description:
 ##     Summary of an object of class dblm.
 ##         - calculate the pearson residuals
 ##         - the deviance residuals
 ##         - the dispersion
 ##        


summary.dbglm <-function(object,dispersion=NULL,...){

   # recover attributs of dbglm object
   y <-object$y
   mu <-object$fitted.values
   family <-object$family
   dev.resids <-family$dev.resids
   df.r <- object$df.residual
   new_weights <-object$weights
   weights <-object$prior.weights
   residuals <-object$residuals
   call<-object$call
   G<-attr(object,"G")
   
   # G vector
   gvec<-diag(G)
   
   # Geometric variability
   gvar=t(weights/sum(weights))%*%as.matrix(gvec)
   
   # pearson residuals
   pears.resid<-(y-mu)/sqrt(family$variance(mu))
   
   # deviance residuals 
   deviance.resid<-sign(y-mu)*sqrt(dev.resids(y,mu,weights))

   # dispersion
    est.disp <- FALSE
    if (is.null(dispersion)){
     dispersion <- if (family$family %in% c("poisson","binomial")) 1
                  else if (df.r > 0) {
                    est.disp <- TRUE
                    if (any(new_weights == 0))
                      warning("observations with zero weight not used for calculating dispersion")
                    sum((new_weights * residuals^2)[new_weights > 0])/df.r
                 }
                 else {
                   est.disp <- TRUE
                   NaN
                 }
   } else{
    if (!is.numeric)
      stop("dispersion must be numeric")
   }

   
  # return a list with the following attributes
  ans<-list(call=call,family=family,deviance=object$deviance,aic=object$aic.model,
            bic = object$bic.model, gcv = object$gcv.model,
            df.residual=df.r,null.deviance=object$null.deviance,
            df.null=object$df.null,iter=object$iter,deviance.resid=deviance.resid,
            pears.resid=pears.resid,dispersion=dispersion,gvar=gvar,gvec=gvec,
            convcrit=object$convcrit)

   class(ans) <- "summary.dbglm"
   return(ans)
}