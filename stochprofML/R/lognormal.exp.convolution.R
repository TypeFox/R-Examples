lognormal.exp.convolution <-
function(z,j.vector,mu.vector,sigma.vector,lambda,logdens) {
# Approximates the density of the sum of a number of independent
# random variables. These are as follows:
#
# - j.vector=(j_1,...j_T), and j_i is the number of summands with distribution i
# - distributions 1 to T-1 are lognormal distributions, distribution T is an exponential distribution
# - mu.vector=(mu_1,...,mu_(T-1)) and sigma.vector=(sigma_1,...,sigma_(T-1)) contain the log-means 
#   and log-standard deviations for the lognormal distributions
# - lambda is the rate of the exponential distribution
# - the density is evaluated at z, which can be multi-dimensional
# - if logdens==T, the log of this density is returned.
#
# Note: 
# If all summands are lognormally distributed, the approximation method by Fenton is used. 
# If all summands are exponentially distributed, their sum follows a gamma distribution.
# If both types of distributions, this function approximates the convolution between the two
# kinds using numerical integration.
#
#
# References:
# [1] L. Fenton (1960). The Sum of Log-Normal Probability Distributions in Scatter
# Transmission Systems. IRE Transactions on Communication Systems, 8: 57-67.


   TY <- length(j.vector)
      
   if (j.vector[TY]==sum(j.vector)) {
      # there are no lognormal types --> easy     
      return(d.sum.of.exp.types(y=z,j=j.vector[TY],lambda=lambda,logdens=logdens))
   }
   else if (j.vector[TY]==0) {
      # there are only lognormal types
      return(d.sum.of.lognormal.types(y=z,j.vector=j.vector[-TY],mu.vector=mu.vector,sigma.vector=sigma.vector,logdens=logdens))
   }


   # otherwise: convolution of lognormal and gamma
   d <- rep(NA,length(z))
   for (i in 1:length(z)) {
      convol <- function(x) {
         # lognormal density
         d1 <- d.sum.of.lognormal.types(y=x,j.vector=j.vector[-TY],mu.vector=mu.vector,sigma.vector=sigma.vector,logdens=F)
         # gamma density
         d2 <- d.sum.of.exp.types(y=z[i]-x,j=j.vector[TY],lambda=lambda,logdens=F)
         
         res <- d1*d2     
         return(res)
      }
      d[i] <- integrate(convol,lower=0,upper=z[i],subdivisions=500,stop.on.error=F)$value
   }
         
   if (logdens) { log(d) } else { d }
}
