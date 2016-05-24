######################################################################
# This file contains utility functions for the generalized Dirichlet
# distribution described in the article by T.-T. Wong et al. (1998),
# Generalized Dirichlet distribution in Bayesian analysis. Applied
# Mathematics and Computation, volume 97, pp 165-181.
#
# This includes:
#  rgd - sample from the generalized Dirichlet distribution
#  Egd - expectation of the generalized Dirichlet distribution
#
# Author: Michael Höhle <hoehle@math.su.se>
# Date:   LaMo Apr 2014. 
######################################################################


######################################################################
# Sample from the generalized dirichlet distribution, i.e.
#  (X_1,...,X_{k+1})' ~ GD(alpha,beta)
# This is the algorithm described by Wong (1998), p. 174.
#
# Parameters:
#  alpha - vector of length k
#  beta  - vector of length k
#
# Note: The alpha and beta vectors are for the first k categories.
# The element in k+1 is automatically given as 1-sum_{i=1}^k X_i.
######################################################################

rgd <- function(n,alpha,beta) {
  #Check that alpha and beta are of the same length.
  if (length(alpha) != length(beta)) {
    stop("alpha and beta not of same length")
  }
  if (!all(alpha>0) | !all(beta>0)) {
    stop("Assumptiom alpha>0 and beta>0 is violated.")
  }
  #Prepare result and sample the first step as in Wong (1998), p.174
  res <- matrix(NA,nrow=n,ncol=length(alpha)+1)
  res[,1] <- rbeta(n,alpha[1],beta[1])
  sum <- res[,1]
  for (j in 2:(length(alpha))) {
    xj <-  rbeta(n, alpha[j], beta[j])
    #Adjust for previous samples
    res[,j] <- xj * (1-sum)
    sum <- sum + res[,j]
  }
  #Last cell is fixed.
  res[,length(alpha)+1] <- 1-sum
  
  return(res) 
}

######################################################################
#Compute analytically the expectation of a GD(alpha,beta) distributed
#variable using the expression of Wong (1998).
#
# Parameters:
#  alpha - vector of alpha parameters of the distribution
#  beta  - vector of beta parameters of the distribution
#
# Returns:
#  Expectation vector of the GD(alpha,betra) distribution
######################################################################
Egd <- function(alpha, beta) {
  mu <- alpha/(alpha+beta)
  mean <- NULL
  for (j in 1:length(mu)) {
    mean[j] <- mu[j] * prod(1-mu[seq_len(j-1)])
  }
  return(c(mean,prod(1-mu)))
}

