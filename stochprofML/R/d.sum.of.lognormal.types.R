d.sum.of.lognormal.types <-
function(y,j.vector,mu.vector,sigma.vector,logdens=F) {
# Density for a sum of independent lognormally distributed variables, 
# using the approximation method by Fenton [1].
#
# - density is evaluated at y, which can be multi-dimensional
#
# - j.vector=(j1,...,jT) is a vector indicating how many of the summands are of
#   which type: 
#         j1 are of type I, ..., jT are of type T.
#   The sum n=j1+...+jT implies how many summands are entering the sum.
#
# - mu.vector=(mu1,...,muT) is of the same length as j.vector. mu_i is the 
#   log-mean for type i.
#
# - sigma.vector is defined analogously as mu.vector. It contains the standard
#   deviations (not variances!) of the according normal distributions.   
# 
# - if logdens==T, the log of this density is returned.
#
#
# References:
# [1] L. Fenton (1960). The Sum of Log-Normal Probability Distributions in Scatter
# Transmission Systems. IRE Transactions on Communication Systems, 8: 57-67.

    
   # check for some obvious errors
   if ((length(j.vector)!=length(mu.vector)) || (length(j.vector)!=length(sigma.vector))) {
      stop("d.sum.of.lognormal.types: j and mu and/or sigma are of different lengths.")
   }
   
   # check whether there are any lognormal summands at all
   if (sum(j.vector)==0) {
      stop("d.sum.of.lognormal.types: sum(j.vector)==0")   
   }   
   
   # these are the "extended" mu and sigma vector, containing 
   # j1 times mu1 (or sigma1, resp), j2 times mu2 (or sigma2, resp), ...
   full.mu.vector <- rep(0,sum(j.vector))
   full.sigma.vector <- full.mu.vector
   
   index <- 0
   for (i in 1:length(j.vector)) {
      if (index<sum(j.vector)) {
         full.mu.vector[index+(1:j.vector[i])] <- mu.vector[i]
         full.sigma.vector[index+(1:j.vector[i])] <- sigma.vector[i]
         index <- index + j.vector[i]
      }
   }
   
   # pass these extended mu and sigma vectors to the actual approximation
   # of the density of sums of lognormals
   d <- d.sum.of.lognormals(y,full.mu.vector,full.sigma.vector)
   if (logdens) { log(d) } else { d }
}
