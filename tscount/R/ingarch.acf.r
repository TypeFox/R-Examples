ingarch.acf <- function(intercept, past_obs=NULL, past_mean=NULL, lag.max=10, type=c("acf", "pacf", "acvf"), plot=TRUE, ...){
#Theoretical autocorrelation function of a Poisson INGARCH(p,q) process
##############################
  ingarch.parametercheck(param=list(intercept=intercept, past_obs=past_obs, past_mean=past_mean, xreg=NULL))  
  type <- match.arg(type)
  p <- length(past_mean)
  q <- length(past_obs)
  if(p<q) past_mean <- c(past_mean, rep(0,q-p)) 
  if(q<p) past_obs <- c(past_obs, rep(0,p-q))
  if(type %in% c("acf", "acvf")){
    acvf <- tacvfARMA(phi=past_obs+past_mean, theta=past_mean, maxLag=lag.max, sigma2=ingarch.mean(intercept=intercept, past_obs=past_obs, past_mean=past_mean)) #Error if not stationary/causal, see help(tacvfARMA)
    result <- switch(type, acf=acvf/acvf[1], acvf=acvf)
    names(result) <- 0:lag.max
  }else{
    result <- ARMAacf(ar=past_obs+past_mean, ma=-past_mean, lag.max=lag.max, pacf=TRUE)
    names(result) <- 1:lag.max
  }  
  if(plot){
    plot(result, type="h", xlab="Lag", ylab=toupper(type), ...)
    abline(h=0)
    return(invisible(result))
  }else{ 
    return(result)
  }
}
