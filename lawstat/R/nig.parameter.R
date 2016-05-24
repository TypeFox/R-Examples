`nig.parameter` <-
function(mean=mean, variance=variance, kurtosis=kurtosis, skewness=skewness)
{
### stop the code if the parameters do not satisfy the constraints. ###
 
 if (3*kurtosis <= 5*skewness^2)
 {
  stop("3*kurtosis-5*skewness^2 must be greater than 0.")
 }

 if (skewness < 0)
 {
  stop("skewness must be nonnegative.")
 }

 if (variance < 0)
 {
  stop("variance must be greater than 0.")
 }

### parameter calculations ###

 alpha<-sqrt(9*(3*kurtosis-4*skewness^2)/(variance*(3*kurtosis-5*skewness^2)^2))
 beta<-sqrt(9*skewness^2/(variance*(3*kurtosis-5*skewness^2)^2))
 delta<-3*sqrt(variance*(3*kurtosis-5*skewness^2))/(3*kurtosis-4*skewness^2)
 mu<-mean-delta*beta/sqrt(alpha^2-beta^2)

### display output ###

 return(list(alpha=alpha,beta=beta,delta=delta,mu=mu))
}

