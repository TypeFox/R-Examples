######################################################################
# Function to calculate a mid-p confidence interval for the
# binomial proportion from one observation of an X \sim Bin(n,pi) 
# distributed random variate.
#
# Params:
#  x - Vector of number of successes in the binomial experiment
#  n - Vector of number of independent trials in the binomial
#          experiment
#  conf.level - The level of confidence to be used in the confidence
#          interval, i.e. conf.level = 1 - alpha
#
# Returns:
#     A data.frame containing the observed proportions and the lower and
#     upper bounds of the confidence interval. The style is similar
#     to the binom.confint function of the binom package
######################################################################

binom.midp <- function(x, n, conf.level = 0.95) {
  #Functions to find root of for the lower and higher bounds of the CI
  #These are helper functions.
  f.low <- function(pi,x,n) {
    1/2*dbinom(x,size=n,prob=pi) + pbinom(x,size=n,prob=pi,lower.tail=FALSE) - (1-conf.level)/2
  }
  f.up <- function(pi,x,n) {
    1/2*dbinom(x,size=n,prob=pi) + pbinom(x-1,size=n,prob=pi) - (1-conf.level)/2
  }

  #Function to calculate the midp confidence interval for one
  #set of (x,n) values
  one.ci <- function(x,n) {
    #One takes pi_low = 0 when x=0 and pi_up=1 when x=n
    pi.low <- 0
    pi.up  <- 1

    #Calculate CI by finding roots of the f funcs
    if (x!=0) {
      pi.low <- uniroot(f.low,interval=c(0,x/n),x=x,n=n)$root
    } 
    if (x!=n) {
      pi.up  <- uniroot(f.up,interval=c(x/n,1),x=x,n=n)$root
    }
    #Done
    return(c(pi.low,pi.up))
  }

  #Handle possible vector arguments
#  if ((length(x) != length(n))) {
    m <- cbind(x = x, n = n)
    x <- m[, "x"]
    n <- m[, "n"]
#  } 

  #Run function
  cis <- t(sapply(1:nrow(m), function(i) one.ci(x=x[i],n=n[i])))
  #Wrap up result
  res <- data.frame(method="midp",x=x,n=n,mean=x/n,lower=cis[,1],upper=cis[,2])

  return(res)
}

#binom.midp(x=0:10,n=10,conf.level=0.95)
