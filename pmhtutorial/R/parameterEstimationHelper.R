##############################################################################
#
# Example of particle Metropolis-Hastings 
#
# Subroutine for particle Metropolis-Hastings
#
# Copyright (C) 2015 Johan Dahlin < johan.dahlin (at) liu.se >
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
##############################################################################

##############################################################################
# Particle Metropolis-Hastings (LGSS model)
##############################################################################

#' Particle Metropolis-Hastings algorithm for a linear Gaussian state space 
#' model
#' @description 
#' Estimates the parameter posterior for \eqn{phi} a linear Gaussian state 
#' space model of the form \eqn{ x_{t} = \phi x_{t-1} + \sigma_v v_t } and 
#' \eqn{ y_t = x_t + \sigma_e e_t }, where \eqn{v_t} and \eqn{e_t} denote 
#' independent standard Gaussian random variables, i.e.\eqn{N(0,1)}.
#' @param y Observations from the model for \eqn{t=1,...,T}.
#' @param initPar The mean of the log-volatility process \eqn{\mu}.
#' @param sigmav The standard deviation of the state process \eqn{\sigma_v}.
#' @param sigmae The standard deviation of the observation process 
#' \eqn{\sigma_e}.
#' @param nPart The number of particles to use in the filter.
#' @param T The number of observations.
#' @param x0 The inital state.
#' @param nIter The number of iterations in the PMH algorithm.
#' @param stepSize The standard deviation of the Gaussian random walk proposal 
#' for \eqn{\phi}.
#'
#' @return
#' The trace of the Markov chain exploring the marginal posterior for 
#' \eqn{\phi}.
#' @references 
#' Dahlin, J. & Schoen, T. B. "Getting started with particle 
#' Metropolis-Hastings for inference in nonlinear dynamical models." 
#' pre-print, arXiv:1511.01707, 2015.
#' @author 
#' Johan Dahlin <johan.dahlin@liu.se>
#' @note 
#' See Section 4.3 in the reference for more details.
#' @keywords 
#' ts
#' @export
#'
#' @example ./examples/pmh
#' @importFrom stats dgamma
#' @importFrom stats dnorm
#' @importFrom stats rnorm
#' @importFrom stats runif

pmh <- function(y, initPar, sigmav, sigmae, nPart, T, x0, nIter, stepSize) {

  #===========================================================
  # Initialise variables
  #===========================================================
  th     <- matrix(0, nrow=nIter, ncol=1)
  thp    <- matrix(0, nrow=nIter, ncol=1)
  ll     <- matrix(0, nrow=nIter, ncol=1)
  llp    <- matrix(0, nrow=nIter, ncol=1)
  accept <- matrix(0, nrow=nIter, ncol=1)
  
  # Set the initial parameter and estimate the initial log-likelihood
  th[1]  <- initPar
  ll[1]  <- sm(y, th[1], sigmav, sigmae, nPart, T, x0)$ll
  
  #=====================================================================
  # Run main loop
  #=====================================================================
  for (kk in 2:nIter) {
    
    # Propose a new parameter
    thp[kk] <- th[kk-1] + stepSize * rnorm(1)
    
    # Estimate the log-likelihood (don't run if unstable system)
    if (abs(thp[kk]) < 1.0) {
      llp[kk] <- sm(y, thp[kk], sigmav, sigmae, nPart, T, x0)$ll
    }
    
    # Compute the acceptance probability
    aprob <- exp(dnorm(thp[kk], log=TRUE) - dnorm(th[kk-1], log=TRUE) + llp[kk] - ll[kk-1])
    
    # Generate uniform random variable in U[0,1]
    u = runif(1)
    
    # Accept / reject step
    # Check if | phi | > 1.0, in that case always reject.
    if ((u < aprob) && ( abs( thp[kk] ) < 1.0 )) {
      # Accept the parameter
      th[kk]     <- thp[kk]
      ll[kk]     <- llp[kk]
      accept[kk] <- 1.0
    } else {
      # Reject the parameter
      th[kk]     <- th[kk-1]
      ll[kk]     <- ll[kk-1]       
      accept[kk] <- 0.0
    }
    
    # Write out progress
    if (kk%%100 == 0) {
      cat(sprintf("#####################################################################\n"))
      cat(sprintf(" Iteration: %d of : %d completed.\n \n", kk, nIter))
      cat(sprintf(" Current state of the Markov chain:       %.4f \n", th[kk] ))
      cat(sprintf(" Proposed next state of the Markov chain: %.4f \n", thp[kk] ))
      cat(sprintf(" Current posterior mean:                  %.4f \n", mean(th[0:kk]) ))
      cat(sprintf(" Current acceptance rate:                 %.4f \n", mean(accept[0:kk]) ))
      cat(sprintf("#####################################################################\n"))
    }
  }
  
  #=====================================================================
  # Return traces of the parameters
  #=====================================================================
  th
}


##############################################################################
# Particle Metropolis-Hastings (SV model)
##############################################################################

#' Particle Metropolis-Hastings algorithm for a stochastic volatility model
#' model
#' @description 
#' Estimates the parameter posterior for \eqn{\theta=\{\mu,\phi,\sigma_v\}} in 
#' a stochastic volatility model of the form \eqn{x_t = \mu + \phi ( x_{t-1} - 
#' \mu ) + \sigma_v v_t} and \eqn{y_t = \exp(x_t/2) e_t}, where \eqn{v_t} and 
#' \eqn{e_t} denote independent standard Gaussian random variables, i.e. 
#' \eqn{N(0,1)}.
#' @param y Observations from the model for \eqn{t=1,...,T}.
#' @param initPar An inital value for the parameters 
#' \eqn{\theta=\{\mu,\phi,\sigma_v\}}.
#' @param nPart The number of particles to use in the filter.
#' @param T The number of observations.
#' @param nIter The number of iterations in the PMH algorithm.
#' @param stepSize The standard deviation of the Gaussian random walk proposal 
#' for \eqn{\theta}.
#'
#' @return
#' The trace of the Markov chain exploring the posterior of \eqn{\theta}.
#' @references 
#' Dahlin, J. & Schoen, T. B. "Getting started with particle 
#' Metropolis-Hastings for inference in nonlinear dynamical models." 
#' pre-print, arXiv:1511.01707, 2015.
#' @author 
#' Johan Dahlin <johan.dahlin@liu.se>
#' @note 
#' See Section 5 in the reference for more details.
#' @keywords 
#' ts
#' @export
#' @example ./examples/pmh_sv
#' @importFrom stats dgamma
#' @importFrom stats dnorm
#' @importFrom stats rnorm
#' @importFrom stats runif

pmh_sv <- function(y, initPar, nPart, T, nIter, stepSize) {

  #===========================================================
  # Initialise variables
  #===========================================================
  xh     <- matrix(0, nrow=nIter, ncol=T+1)
  xhp    <- matrix(0, nrow=nIter, ncol=T+1)
  th     <- matrix(0, nrow=nIter, ncol=3)
  thp    <- matrix(0, nrow=nIter, ncol=3)
  ll     <- matrix(0, nrow=nIter, ncol=1)
  llp    <- matrix(0, nrow=nIter, ncol=1)
  accept <- matrix(0, nrow=nIter, ncol=1)
  
  # Set the initial parameter and estimate the initial log-likelihood
  th[1,]  <- initPar
  res     <- sm_sv(y, th[1,1], th[1,2], th[1,3], nPart, T)
  ll[1]   <- res$ll
  xh[1,]  <- res$xh
  
  #=====================================================================
  # Run main loop
  #=====================================================================
  for (kk in 2:nIter) {
    
    # Propose a new parameter
    thp[kk,] <- mvtnorm::rmvnorm(1, mean=th[kk-1,], sigma=stepSize)
    
    # Estimate the log-likelihood (don't run if unstable system)
    if ((abs(thp[kk,2]) < 1.0 ) && (thp[kk,3] > 0.0 )) {
      res      <- sm_sv(y, thp[kk,1], thp[kk,2], thp[kk,3], nPart, T)
      llp[kk]  <- res$ll
      xhp[kk,] <- res$xh
    }
    
    # Compute difference in the log-priors
    dpmu  <- dnorm(thp[kk,1], 0, 1, log=TRUE) - dnorm(th[kk-1,1], 0, 1, log=TRUE)
    dpphi <- dnorm(thp[kk,2], 0.95, 0.05, log=TRUE) - dnorm(th[kk-1,2], 0.95, 0.05, log=TRUE)
    dpsiv <- dgamma(thp[kk,3], 2, 10, log=TRUE) - dgamma(th[kk-1,3], 2, 10, log=TRUE)
    
    # Compute the acceptance probability
    aprob <- exp(dpmu + dpphi + dpsiv + llp[kk] - ll[kk-1])
    
    # Generate uniform random variable in U[0,1]
    u <- runif(1)
    
    # Accept / reject step
    # Check if | phi | > 1.0, in that case always reject.
    if ((u < aprob) && (abs(thp[kk,2]) < 1.0)) {
      # Accept the parameter
      th[kk,]     <- thp[kk,]
      ll[kk]      <- llp[kk]
      xh[kk,]     <- xhp[kk,]
      accept[kk]  <- 1.0
    } else {
      # Reject the parameter
      th[kk,]     <- th[kk-1,]
      ll[kk]      <- ll[kk-1]      
      xh[kk,]     <- xh[kk-1,]
      accept[kk]  <- 0.0
    }
    
    # Write out progress
    if (kk%%100 == 0) {
      cat(sprintf("#####################################################################\n"))
      cat(sprintf(" Iteration: %d of : %d completed.\n \n", kk, nIter));
      cat(sprintf(" Current state of the Markov chain:       %.4f %.4f %.4f \n", th[kk,1], th[kk,2], th[kk,3] ))
      cat(sprintf(" Proposed next state of the Markov chain: %.4f %.4f %.4f \n", thp[kk,1], thp[kk,2], thp[kk,3] ))
      cat(sprintf(" Current posterior mean:                  %.4f %.4f %.4f \n", mean(th[0:kk,1]), mean(th[0:kk,2]), mean(th[0:kk,3]) ))
      cat(sprintf(" Current acceptance rate:                 %.4f \n", mean(accept[0:kk]) ))
      cat(sprintf("#####################################################################\n"))
    }
  }
  
  #=====================================================================
  # Return traces of the parameters
  #=====================================================================
  list(thhat=th, xhat=xh)
}

##############################################################################
# Particle Metropolis-Hastings (SV model)
##############################################################################

#' Particle Metropolis-Hastings algorithm for a stochastic volatility model
#' model
#' @description 
#' Estimates the parameter posterior for \eqn{\theta=\{\mu,\phi,\sigma_v\}} in 
#' a stochastic volatility model of the form \eqn{x_t = \mu + \phi ( x_{t-1} - 
#' \mu ) + \sigma_v v_t} and \eqn{y_t = \exp(x_t/2) e_t}, where \eqn{v_t} and 
#' \eqn{e_t} denote independent standard Gaussian random variables, i.e. 
#' \eqn{N(0,1)}. In this version of the PMH, we reparameterise the model and 
#' run the Markov chain on the  parameters \eqn{\vartheta=\{\mu,\psi,
#' \varsigma\}}, where \eqn{\phi=\tanh(\psi)} and \eqn{sigma_v=\exp(\varsigma)}.
#' @param y Observations from the model for \eqn{t=1,...,T}.
#' @param initPar An inital value for the parameters 
#' \eqn{\theta=\{\mu,\phi,\sigma_v\}}.
#' @param nPart The number of particles to use in the filter.
#' @param T The number of observations.
#' @param nIter The number of iterations in the PMH algorithm.
#' @param stepSize The standard deviation of the Gaussian random walk proposal 
#' for \eqn{\theta}.
#'
#' @return
#' The trace of the Markov chain exploring the posterior of \eqn{\theta}.
#' @references 
#' Dahlin, J. & Schoen, T. B. "Getting started with particle 
#' Metropolis-Hastings for inference in nonlinear dynamical models." 
#' pre-print, arXiv:1511.01707, 2015.
#' @author 
#' Johan Dahlin <johan.dahlin@liu.se>
#' @note 
#' See Section 5 in the reference for more details.
#' @keywords 
#' ts
#' @export
#' @example ./examples/pmh_sv_reparam
#' @importFrom stats dgamma
#' @importFrom stats dnorm
#' @importFrom stats rnorm
#' @importFrom stats runif

pmh_sv_reparameterised <- function(y, initPar, nPart, T, nIter, stepSize){

  #===========================================================
  # Initialise variables
  #===========================================================
  xh     <- matrix(0, nrow=nIter, ncol=T+1)
  xhp    <- matrix(0, nrow=nIter, ncol=T+1)
  th     <- matrix(0, nrow=nIter, ncol=3)
  thp    <- matrix(0, nrow=nIter, ncol=3)
  tht    <- matrix(0, nrow=nIter, ncol=3)
  thpt   <- matrix(0, nrow=nIter, ncol=3)  
  ll     <- matrix(0, nrow=nIter, ncol=1)
  llp    <- matrix(0, nrow=nIter, ncol=1)
  accept <- matrix(0, nrow=nIter, ncol=1)
  
  # Set the initial parameter and estimate the initial log-likelihood
  tht[1,] <- initPar
  res     <- sm_sv(y, tht[1,1], tht[1,2], tht[1,3], nPart, T)
  th[1,]  <- c(tht[1,1], atanh( tht[1,2]), log(tht[1,3]))
  ll[1]   <- res$ll
  xh[1,]  <- res$xh
  
  #=====================================================================
  # Run main loop
  #=====================================================================
  for (kk in 2:nIter) {
    
    # Propose a new parameter
    thp[kk,] <- mvtnorm::rmvnorm(1, mean=th[kk-1,], sigma=stepSize);
    
    # Run the particle filter
    thpt[kk,]   <- c(thp[kk,1], tanh(thp[kk,2]), exp(thp[kk,3]))
    res         <- sm_sv(y, thpt[ kk,1], thpt[ kk,2], thpt[ kk,3], nPart, T)
    xhp[kk,]    <- res$xh
    llp[ kk ]   <- res$ll
    
    # Compute the acceptance probability
    logPrior1 <- dnorm(thpt[kk,1], log=TRUE) - dnorm(tht[kk-1,1], log=TRUE)
    logPrior2 <- dnorm(thpt[kk,2], 0.95, 0.05, log=TRUE) - dnorm(tht[kk-1,2], 0.95, 0.05, log=TRUE)
    logPrior3 <- dgamma(thpt[kk,3], 3, 10, log=TRUE) - dgamma(tht[ kk-1, 3 ], 3, 10, log=TRUE)
    
    logJacob  <- log(abs(1 - thpt[kk,2 ]^2)) - log(abs(1 - tht[kk-1,2]^2)) + log(abs(thpt[kk,3])) - log(abs(tht[kk-1,3]))
    aprob     <- exp(logPrior1 + logPrior2 + logPrior3 + llp[kk] - ll[kk-1] + logJacob)
    
    # Generate uniform random variable in U[0,1]
    u <- runif(1)
    
    # Accept / reject step
    if (u < aprob) {
      # Accept the parameter
      th[kk,]     <- thp[kk,]
      tht[kk,]    <- thpt[kk,]
      ll[kk]      <- llp[kk]
      xh[kk,]     <- xhp[kk,]
      accept[kk]  <- 1.0
    } else {
      # Reject the parameter
      th[kk,]     <- th[kk-1,]
      tht[kk,]    <- tht[kk-1,]
      ll[kk]      <- ll[kk-1]      
      xh[kk,]     <- xh[kk-1,]
      accept[kk]  <- 0.0
    }
    
    # Write out progress
    if (kk%%100 == 0) {
      cat(sprintf("#####################################################################\n"))
      cat(sprintf(" Iteration: %d of : %d completed.\n \n", kk, nIter))
      cat(sprintf(" Current state of the Markov chain:       %.4f %.4f %.4f \n", tht[kk,1], tht[kk,2], tht[kk,3] ))
      cat(sprintf(" Proposed next state of the Markov chain: %.4f %.4f %.4f \n", thpt[kk,1], thpt[kk,2], thpt[kk,3] ))
      cat(sprintf(" Current posterior mean:                  %.4f %.4f %.4f \n", mean(tht[0:kk,1]), mean(tht[0:kk,2]), mean(tht[0:kk,3]) ))
      cat(sprintf(" Current acceptance rate:                 %.4f \n", mean(accept[0:kk]) ))
      cat(sprintf("#####################################################################\n"));
    }
  }
  
  #=====================================================================
  # Return traces of the parameters
  #=====================================================================
  list(thhat=tht, xhat=xh, thhattansformed=th)
}
##############################################################################
# End of file
##############################################################################
