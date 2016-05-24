##########################################################################
## simple instructional models using Monte Carlo simulation
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Copyright (C) 2003-2007 Andrew D. Martin and Kevin M. Quinn
## Copyright (C) 2007-present Andrew D. Martin, Kevin M. Quinn,
##    and Jong Hee Park
##########################################################################

## Monte Carlo simulation from the likelihood of a 
## binomial distribution with a Beta(alpha, beta) prior
## ADM 1/25/2006


MCbinomialbeta <- function(y, n, alpha=1, beta=1, mc=1000, ...) {
    
   # check data
   if(y < 0) {
      cat("Error: Number of successes negative.\n")
      stop("Please respecify and call function again.")
    }
   if(n < 0) {
      cat("Error: Number of trials negative.\n")
      stop("Please respecify and call function again.")  
   }
   if(y > n) {
      cat("Error: Number of successes greater than number of trials.\n")
      stop("Please respecify and call function again.")   
   }
   
   # check other parameters
   check.beta.prior(alpha, beta)
   check.mc.parameter(mc)
   
   # draw sample and return
   output <- mcmc(matrix(rbeta(mc, alpha+y, beta+n-y),mc,1))
   varnames(output) <- as.list("pi")
   attr(output,"title") <- "MCbinomialbeta Posterior Sample"
      
   return(output)
}

## Monte Carlo simulation from the likelihood of a 
## Poisson distribution with a Gamma(alpha, beta) prior
## ADM 1/25/2006
MCpoissongamma <- function(y, alpha, beta, mc=1000, ...) {
    
   # check data
   if(any(y < 0)) {
      cat("Error: Some counts negative in y.\n")
      stop("Please respecify and call function again.")
    }
   n <- length(y)
   
   # check other parameters
   check.gamma.prior(alpha, beta)
   check.mc.parameter(mc)
   
   # draw sample and return
   output <- mcmc(matrix(rgamma(mc, alpha+sum(y), beta+n),mc,1))
   varnames(output) <- as.list("lambda")
   attr(output,"title") <- "MCpoissongamma Posterior Sample"
      
   return(output)
}

## Monte Carlo simulation from the likelihood of a 
## Normal distribution with a Normal(mu0, tau20) prior
## the variance sigma2 is known
## ADM 1/26/2006
MCnormalnormal <- function(y, sigma2, mu0, tau20, mc=1000, ...) {
   
   n <- length(y)
   if(sigma2 <= 0) {
      cat("Error: Known variance sigma2 is less than or equal to zero.\n")
      stop("Please respecify and call function again.")
    }
   
   # check other parameters
   check.normal.prior(mu0, tau20)
   check.mc.parameter(mc)
   
   # draw sample and return
   mu1 = (1/tau20 * mu0 + n/sigma2 * mean(y)) / (1/tau20 + n/sigma2)
   tau21 = 1 / (1/tau20 + n/sigma2)
   output <- mcmc(matrix(rnorm(mc, mu1, sqrt(tau21)),mc,1))
   varnames(output) <- as.list("mu")
   attr(output,"title") <- "MCnormalnormal Posterior Sample"

   return(output)
}

## Monte Carlo simulation from the likelihood of a 
## multinomal distribution with a Dirichlet(alpha) prior
MCmultinomdirichlet <- function(y, alpha0, mc=1000, ...) {
   
   # check data
   d <- length(y)
   if(any(y < 0)) {
      cat("Error: Some counts negative in y.\n")
      stop("Please respecify and call function again.")
    }
   
   # check alpha
   if(length(alpha0) != d) {
      cat("Error: Dimension of alpha and y do not match.\n")
      stop("Please respecify and call function again.")    
   }
   if(any(alpha0 <= 0)) {
      cat("Error: At least one alpha in Dirichlet prior less than or equal to zero.\n")
      stop("Please respecify and call function again.")   
   }
   
   # draw sample and return
   output <- mcmc(rdirichlet(mc,y + alpha0))
   varnames(output) <- paste("pi.", 1:d, sep="")
   attr(output,"title") <- "MCmultinomdirichlet Posterior Sample"

   return(output)
}
