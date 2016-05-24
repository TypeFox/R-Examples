interv_test <- function(...) UseMethod("interv_test")

interv_test.tsglm <- function(fit, tau, delta, external, info=c("score"), est_interv=FALSE, ...){
#Test on one or several interventions of known types at known points in time
##############################

  #Check and modify argument:  
  tsglm.check(fit)
  info <- match.arg(info)
  if(info=="hessian" && fit$link=="link") stop("For a model with logarithmic link argument 'info' needs to be set to \"score\"")
  r <- length(tau) #number of additional parameters
  if(missing(external) || length(external)==0) external <- rep(FALSE, r) else external <- as.logical(external) #the default value for external is FALSE (i.e. an internal intervention effect)
  if(length(external)==1) external <-  rep(external, r) else external <- as.logical(external) #if only one value for external is provided, this is used for all interventions
  
  #Add information about intervention effects:
    param_H0_extended <- c(fit$coefficients, numeric(r))
    model_extended <- fit$model
    xreg_extended <- fit$xreg
    covariate <- interv_covariate(n=length(fit$ts), tau=tau, delta=delta)
    xreg_extended <- cbind(fit$xreg, covariate)
    colnames(xreg_extended) <- c(colnames(fit$xreg), colnames(covariate)) 
    model_extended$external <- c(fit$model$external, external)
    loglik <- tsglm.loglik(link=fit$link, paramvec=param_H0_extended, model=model_extended, ts=fit$ts, xreg=xreg_extended, score=TRUE, info=info)   
  infomat_corrected <- apply((1/loglik$kappa + fit$sigmasq)*loglik$outerscoreprod, c(2,3), sum)
  test_statistic <- scoretest(Score=loglik$score, G=loglik$info, G1=infomat_corrected, r=r, stopOnError=TRUE)$test_statistic    
  p_value <- 1-pchisq(test_statistic, df=r)
  result <- list(
    test_statistic=test_statistic,
    df=r,
    p_value=p_value,
    fit_H0=fit,
    model_interv=model_extended,
    xreg_interv=xreg_extended
  )
  if(est_interv){ #ML estimation for the model with intervention
      fit_interv <- tsglm(ts=fit$ts, model=model_extended, xreg=xreg_extended, link=fit$link, distr=fit$distr, ...)
      result <- c(result, list(fit_interv=fit_interv))  
  }
  class(result) <- "interv_test"
  return(result)
}
