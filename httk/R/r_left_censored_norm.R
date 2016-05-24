# This function returns draws from a normal distribution with a lower censoring
# limit of lod (limit of detection):
r.left.censored.norm <- function(n,mean=0,sd=1,lod=0.5)
{
#  erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1
#  lod.CDF <- 0.5*(1+erf((lod-mean)/(2*sd^2)^(1/2)))
  censored <- runif(n,0,1) < ptnorm(lod,mean=mean,sd=sd,lower=0)
  out <- rep(0,length=n)
  out[censored] <- runif(sum(censored),0,lod)
  out[!censored] <- rtnorm(sum(!censored),mean=mean,sd=sd,lower=lod)
  return(out)
}