ingarch.parametercheck <- function(param){
  stopifnot(     
    param$intercept>0,
    param$past_obs>=0,
    param$past_mean>=0,
    param$xreg>=0
  )
  sum_param_past <- sum(param$past_obs)+sum(param$past_mean)
  if(sum_param_past>=1) stop(paste("Parameters are outside the stationary region, sum of parameters for regression\non past observations and on past conditional means is", sum_param_past, "> 1"))
}
 