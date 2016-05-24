d.sum.of.lognormals <-
function(y,mu.vector,sigma.vector) {
# Approximates the density of the sum of a number of independent
# lognormally distributed random variables. Since this is not possible in 
# analytically closed form, the approximation by Fenton [1] is used.
#
# The length of mu.vector (and sigma.vector) is the number of random
# variables entering the sum. The i.th element of mu.vector
# (and sigma.vector) is the log-mean (and log-standard deviation) of the i.th
# random variable. The density is evaluated at y.
#
# References:
# [1] L. Fenton (1960). The Sum of Log-Normal Probability Distributions in Scatter
# Transmission Systems. IRE Transactions on Communication Systems, 8: 57-67.


   # calculate the mean and variance of the sum of lognormal distributions
   E <- sum(exp(mu.vector + 0.5*sigma.vector^2))
   log.Vi.part1 <- 2*mu.vector
   log.Vi.part2 <- sigma.vector^2
   log.Vi.part3 <- log(exp(sigma.vector^2)-1)
   
   if (sum(abs(log.Vi.part3)==Inf)>0) {
      # the second or third part will be dominating
      V <- sum(exp(log.Vi.part3))
   }
   else {
      V <- sum(exp(log.Vi.part1+log.Vi.part2+log.Vi.part3))
   }
       
      
   # check whether these parameters are "proper":
   #
   # if the variance is zero, the density is a point mass around the expectation
   if (round(V,4)==0) {  
      k <- which(y==E)
      d <- rep(0,length(y))
      d[k] <- 10^7 # use finite value instead of Inf because optim does not like Inf
      return(d)
   }
   # if the variance is very large or the expectation is very large or zero,
   # the density is numerically equal to zero
   else if (((V>10^10) || (E>10^10)) || (E==0)) {
      return(rep(0,length(y)))
   }
   
   # calculate the mean and variance of the according normal distribution   
   m <- log(E)-0.5*log(V/E^2+1)  
   s2 <- log(V/E^2+1)
      
   # if the variance is infinite, the density is equal to zero   
   if (s2==Inf) {
      return(rep(0,length(y)))
   }
   # if the variance is zero, the density is a point mass around the expectation
   if (round(sqrt(s2),4)==0) {  
      k <- which(y==E)
      d <- rep(0,length(y))
      d[k] <- 10^7 # use finite value instead of Inf because optim does not like Inf
      return(d)
   }
      
   # The sum of lognormal variables is approximately lognormally distributed 
   # with log-mean m and log-standard deviation sqrt(s2).
   return(dlnorm(y,m,sqrt(s2)))
}
