######################################################################
# Compute LRT based intervals for the binomial response
#  W \sim Bin(m, theta), where theta = 1 - (1-pi)^k.
#  Hence pi = g(theta) = 1 - (1-pi)^(1/k)
# One then knows that the borders for pi are just transformations
# of the borders of theta using the above g(theta) function.
#
# Params:
#   ...
######################################################################

poolbinom.lrt <- function(x, k, n, conf.level = 0.95, bayes=FALSE, conf.adj=FALSE) {

  #Compute LRT intervals for binomial response with parameter theta
  cis.theta <- binom.lrt(x=x, n=n, conf.level=conf.level, bayes=bayes, conf.adj=conf.adj, plot=FALSE)

  #Transformation function between observed Bin(m,theta) and pi of the individual Bernoulli 
  g <- function(pi, x, k)  return(1 - (1- pi)^(1/k))
  
  #Transform borders into ci borders for pi according to formula pi = 1 - (1- theta)^(1/k)
  cis <- g(cis.theta[,c("lower","upper")],x=x,k=k)

  #Wrap up result
  res <- data.frame(method="logit",x=x,k=k, n=n, mle=g(cis.theta$mean,x=x,k=k),lower=cis[,1],upper=cis[,2])
  return(res)
}

#debug(poolbinom.lrt)
#debug(ciss.binom)
#poolbinom.

