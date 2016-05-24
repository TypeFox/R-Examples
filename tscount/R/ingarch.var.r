ingarch.var <- function(intercept, past_obs=NULL, past_mean=NULL){
#Theoretical marginal variance of a Poisson INGARCH(p,q) process
##############################
  ingarch.parametercheck(param=list(intercept=intercept, past_obs=past_obs, past_mean=past_mean, xreg=NULL))
  result <- ingarch.acf(intercept=intercept, past_obs=past_obs, past_mean=past_mean, lag.max=0, type="acvf", plot=FALSE)[[1]]
  return(result)
}
