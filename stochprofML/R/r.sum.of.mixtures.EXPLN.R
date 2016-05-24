r.sum.of.mixtures.EXPLN <-
function(k,n,p.vector,mu.vector,sigma.vector,lambda) {
# Draws k i.i.d. random variables from the following distribution:
# Each random variable is the sum of another n independent random variables.
# These are from a mixture of zero, one or more lognormal distributions and
# one exponential distribution. More specifically, with probability p_i (i=1,...,T-1), 
# a summand is of type i. In that case, it is lognormally 
# distributed with log-mean mu_i and log-standard deviation sigma_i.
# With probability p_T, it is exponentially distributed with rate lambda.
#
# Parameters:
#
# - k is the number of i.i.d. random variables returned by this function
#   (in the considered application: number of tissue samples)
# - n is the number of summands entering each of the k random variables
#   (in the considered application: number of cells per tissue sample)
# - p.vector=(p1,p2,..,pT) is the probability for each type
# - mu.vector=(mu1,mu2,...,mu{T-1]) is the log-mean for each lognormal type
# - sigma.vector=(sigma1,...,sigma{T-1}) is the log-standard deviation for each lognormal type     
# - lambda is the rate for the exponential type
#
# The lengths mu.vector and sigma.vector have to be identical. 
# p.vector has to have one component more.
# lambda has to be a scalar.
# These lengths automatically determine the number of different types.

   # check for some obvious errors
   if (round(sum(p.vector),4)!=1) {
      stop("r.sum.of.mixtures: Sum of p's does not equal one.")
   }
   if (sum(p.vector>1)+sum(p.vector<0)>0) {
      stop("r.sum.of.mixtures: There are p's which are not in [0,1]")
   }
   if ((length(p.vector)!=length(mu.vector)+1) || (length(p.vector)!=length(sigma.vector)+1)) {
      stop("r.sum.of.mixtures: p and mu and/or sigma are of contradicting lengths.")
   }
   if (length(lambda)!=1) {
      stop("r.sum.of.mixtures: lambda is not a scalar.")
   }

   TY <- length(p.vector)

   # draw the number of summands of type i (i.th column) in each sample j (j.th row)
   N.matrix <- t(rmultinom(n=k, size=n, prob=p.vector))
   
   # draw the summands and sum up
   random <- rep(NA,k)
   for (j in 1:k) {   
      # j.th sample
      N.vector <- N.matrix[j,]
      # expand mu.vector and sigma.vector
      full.mu.vector <- rep(0,n-N.vector[TY])
      full.sigma.vector <- full.mu.vector
      index <- 0
      for (i in 1:length(p.vector)) {
         if (index<n-N.vector[TY]) {
            full.mu.vector[index+(1:N.vector[i])] <- mu.vector[i]
            full.sigma.vector[index+(1:N.vector[i])] <- sigma.vector[i]
            index <- index + N.vector[i]
         }
      }
      # vector of lognormal summands
      Y.ln <- rlnorm(n-N.vector[TY],full.mu.vector,full.sigma.vector)
      # vector of exponential summands
      Y.exp <- rexp(N.vector[TY],lambda)
      # sum up
      random[j] <- sum(Y.ln) + sum(Y.exp)
   }   
   
   return(random)
}
