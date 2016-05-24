########################################################################
#This function calculates the lower/upper CI limit given the estimation, variance and significance level
###########################################
CI.normal <-
function(est,var0,alpha=0.05)
{
  z0 <- abs(qnorm(alpha/2))
  lower <- est-z0*sqrt(var0)
  upper <- est+z0*sqrt(var0)  
  return(list(est=est,var0=var0,lower=lower,upper=upper))
}

