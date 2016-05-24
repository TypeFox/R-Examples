summary.tsglm <- function(object, B, parallel=FALSE, ...){
  cl <- object$call
  if(length(coef(object))){
    ll <- logLik(object) #log-likelihood of the fitted model
    k <- attr(ll, "df") #number of parameters of the fitted model
    n <- attr(ll, "nobs") #number of observations used for model fit
    infer <- NULL
    try(infer <- se(object, B=B, parallel=parallel))  
    if(is.null(object$info.matrix_corrected) | is.null(infer)){
      coefs <- data.frame(Estimate=c(coef(object), object$distrcoefs))
      se_info <- NULL
    }else{
      coefs <- data.frame(
        infer$est,
        infer$se
      )
      names(coefs) <- c("Estimate", "Std. Error")
      se_info <- list(se.type=infer$type)
      if(infer$type=="bootstrap") se_info <- c(se_info, list(se.bootstrapsamples=infer$B))
    }
    result <- c(
      list(
        call=cl,
        link=object$link,
        distr=object$distr,
        residuals=residuals(object, type="response"),
        coefficients=coefs,
        number.coef=nrow(coefs)
      ),
      se_info,            
      list(
        logLik=logLik(object),             
        AIC=AIC(object), #Akaike's Information Criterion
        BIC=BIC(object), #Bayesian Information Criterion
        pearson.resid=residuals(object, type="pearson") #Pearson's residuals
      )
    )  
  }else{ 
    result <- list(call=cl, init=object$init)
  }    
  class(result) <- "summary.tsglm"
  return(result)
}
