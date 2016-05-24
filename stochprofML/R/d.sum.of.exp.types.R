d.sum.of.exp.types <-
function(y,j,lambda,logdens=F) {
# Density for a sum of independent exponentially distributed variables. 
#
# - density is evaluated at y, which can be multi-dimensional
#
# - j is the number of summands.
#
# - lambda is the rate of the exponential distibution.
# 
# - if logdens==T, the log of this density is returned.
    
   # check for some obvious errors
   if (lambda < 0) {
      print(lambda)
      stop("d.sum.of.exp.types: lambda is non-positive.")
   }  
   # check whether there are any lognormal summands at all
   if (j==0) {
      stop("d.sum.of.exp.types: j==0")
   }
   
   # density

   if (lambda==Inf) {
      # since 0 < j < Inf, the density is then equal to Dirac(y)
      d <- rep(0,length(y))
      d[y==0] <- 10^7 # use finite value instead of Inf because optim does not like Inf
      if (logdens) { d <- log(d) }
   }
   else {
      d <- dgamma(y,shape=j,rate=lambda,log=logdens)
   }
   return(d)
}
