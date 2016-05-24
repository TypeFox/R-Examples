######################################################################
# Function to calculate confidence interval based on equation (3.1)
# in Liu and Bailey together with the stretching proposed
# on the top of page 3 in their paper.
#
# Params:
#  x - Vector of number of successes in the binomial experiment
#  n - Vector of number of independent trials in the binomial
#          experiment
#  lambda - shrinkage factor for the proportion estimate
#  d - width of the confidence interval
#
#  conf.level - The level of confidence to be used in the confidence
#          interval, i.e. conf.level = 1 - alpha
#
# Returns:
#     A data.frame containing the observed proportions and the lower and
#     upper bounds of the confidence interval. The style is similar
#     to the binom.confint function of the binom package
######################################################################

binom.liubailey <- function(x, n, d, lambda=0, conf.level = 0.95) {
  #Handle possible vector arguments
  if ((length(x) != length(n))) {
    m <- cbind(x = x, n = n)
    x <- m[, "x"]
    n <- m[, "n"]
  }

  #Compute confidence interval
  pi.hat <- x/n
  z <- -qnorm( (1-conf.level)/2)
  #Estimator with shrinkage
  Cnpi.hat <- pi.hat + lambda*z^2*(0.5-pi.hat)/(n+z^2)
  pi.low <- Cnpi.hat - d
  pi.up  <- Cnpi.hat + d

  #Stretch (in the general form)
#  pi.low.star <- pmax(0, pmin(1-2*d, (pi.low+pi.up)/2-d))
#  pi.up.star <-  pmin(1, pmax(2*d, (pi.low+pi.up))/2+d)
  #Stretch for the (3.1) form. 
  pi.low.star <- pmax(0, pmin(1-2*d, pi.low))
  pi.up.star <-  pmin(1, pmax(2*d, pi.up))
  
  #The result in stretched form
  cis <- cbind(pi.low.star,pi.up.star)

  #Wrap up result in a binom like syntax
  res <- data.frame(method="liubailey",x=x,n=n,mean=x/n,lower=cis[,1],upper=cis[,2])

  return(res)
}

