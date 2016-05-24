##############################################################################
#
# Example of particle filtering 
#
# Subroutine for data generation and particle filtering
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
# Generate data for LGSS model
##############################################################################

#' Generates data from a linear Gaussian state space model
#' @description 
#' Generates data from a specific linear Gaussian state space model of the form 
#' \eqn{ x_{t} = \phi x_{t-1} + \sigma_v v_t } and \eqn{ y_t = x_t + 
#' \sigma_e e_t }, where \eqn{v_t} and \eqn{e_t} denote independent standard 
#' Gaussian random variables, i.e. \eqn{N(0,1)}.
#' @param phi The parameter \eqn{\phi} that scales the current state in the 
#' state dynamics. It is restricted to [-1,1] to obtain a stable model.
#' @param sigmav The standard deviation of the state process noise. Must be 
#' positive.
#' @param sigmae The standard deviation of the observation process noise. Must 
#' be positive.
#' @param T The number of time points to simulate.
#' @param x0 The initial state.
#'
#' @return
#' The function returns a list with the elements: 
#' \itemize{
#' \item{x: The latent state for \eqn{t=0,...,T}.}
#' \item{y: The observation for \eqn{t=0,...,T}.}
#' }
#' @references 
#' Dahlin, J. & Schoen, T. B. "Getting started with particle 
#' Metropolis-Hastings for inference in nonlinear dynamical models." 
#' pre-print, arXiv:1511.01707, 2015.
#' @author 
#' Johan Dahlin <johan.dahlin@liu.se>
#' @keywords 
#' datagen
#' @export
#' @importFrom stats rnorm

generateData <- function(phi, sigmav, sigmae, T, x0) {

  # Pre-allocate vectors for log-volatility/state (x) 
  # and log-returns/observations (y)
  x    <- matrix(0, nrow=T+1, ncol=1)
  y    <- matrix(0, nrow=T+1, ncol=1)
  
  # Set the initial state
  x[1] <- x0
  y[1] <- NA
  
  # Simulate the system for each time step
  for (tt in 2:(T+1)) {
    x[tt] <- phi  * x[tt-1] + sigmav * rnorm(1)
    y[tt] <-        x[tt]   + sigmae * rnorm(1)
  }
  
  list(x=x, y=y)
}


##############################################################################
# Fully-adapted particle filter (LGSS)
##############################################################################

#' Fully-adapted particle filter for state estimate in a linear Gaussian state 
#' space model
#' @description 
#' Estimates the filtered state and the log-likelihood for a linear Gaussian 
#' state space model of the form \eqn{ x_{t} = \phi x_{t-1} + \sigma_v v_t } 
#' and \eqn{ y_t = x_t + \sigma_e e_t }, where \eqn{v_t} and \eqn{e_t} denote 
#' independent standard Gaussian random variables, i.e.\eqn{N(0,1)}.
#' @param y Observations from the model for \eqn{t=1,...,T}.
#' @param phi The persistence of the state process \eqn{\phi}.
#' @param sigmav The standard deviation of the state process \eqn{\sigma_v}.
#' @param sigmae The standard deviation of the observation process \eqn{\sigma_e}.
#' @param nPart The number of particles to use in the filter.
#' @param T The number of observations.
#' @param x0 The initial state.
#'
#' @return
#' The function returns a list with the elements:
#' \itemize{
#' \item{xh: The estimate of the filtered state at time \eqn{t=1,...,T}.}
#' \item{ll: The estimate of the log-likelihood.}
#' \item{p: The particle system at each time point.}
#' \item{w: The particle weights at each time point.}
#' }
#' @references 
#' Dahlin, J. & Schoen, T. B. "Getting started with particle 
#' Metropolis-Hastings for inference in nonlinear dynamical models." 
#' pre-print, arXiv:1511.01707, 2015.
#' @author 
#' Johan Dahlin <johan.dahlin@liu.se>
#' @note 
#' See Section 4 in the reference for more details.
#' @keywords 
#' ts
#' @export
#' @example ./examples/sm
#' @importFrom stats dnorm
#' @importFrom stats rnorm

sm <- function(y, phi, sigmav, sigmae, nPart, T, x0) {

  #===========================================================
  # Initialise variables
  #===========================================================
  xhatf <- matrix(x0, nrow=T, ncol=1)
  p     <- matrix(x0, nrow=nPart, ncol=T+1)
  w     <- matrix(1/nPart, nrow=nPart, ncol=T+1)
  v     <- matrix(1, nrow=nPart, ncol=T+1)
  ll    <- 0

  #===========================================================
  # Run main loop
  #===========================================================
  for (tt in 2:T) {
    
    #=========================================================
    # Resample ( multinomial )
    #=========================================================
    nIdx   <- sample(1:nPart, nPart, replace=TRUE, prob = w[,tt-1])
    
    #=========================================================
    # Propagate
    #=========================================================
    Delta  <- ( sigmav^(-2) + sigmae^(-2) )^(-1)
    mup    <- sigmae^(-2) * y[tt] + sigmav^(-2) * phi * p[nIdx,tt-1]
    p[,tt] <- Delta * mup + rnorm(nPart, 0, sqrt(Delta))
    
    #=========================================================
    # Compute weights
    #=========================================================
    v[,tt] <- dnorm(y[tt+1], phi * p[,tt], sqrt( sigmae^2 + sigmav^2), log=TRUE)
    
    # Rescale log-weights and recover weight
    vmax   <- max(v[,tt])
    v[,tt] <- exp(v[,tt] - vmax)
    
    # Normalize the weights
    w[,tt] <- v[,tt] / sum(v[,tt])
        
    # Estimate the state
    xhatf[tt] <- mean(p[,tt])
    
    # Estimate the log-likelihood
    ll        <- ll + vmax + log(sum(v[,tt])) - log(nPart)
    
  }
  
  #===========================================================
  # Return state estimate and log-likelihood estimate
  #===========================================================
  list(xh = xhatf, ll=ll, p=p, w=w)
  
}

###################################################################################
# Kalman filter (LGSS)
###################################################################################

#' Kalman filter for state estimate in a linear Gaussian state space model
#' @description 
#' Estimates the filtered state and the log-likelihood for a linear Gaussian 
#' state space model of the form \eqn{ x_{t} = \phi x_{t-1} + \sigma_v v_t } 
#' and \eqn{ y_t = x_t + \sigma_e e_t }, where \eqn{v_t} and \eqn{e_t} denote 
#' independent standard Gaussian random variables, i.e.\eqn{N(0,1)}.
#' @param y Observations from the model for \eqn{t=1,...,T}.
#' @param phi The persistence of the state process \eqn{\phi}.
#' @param sigmav The standard deviation of the state process \eqn{\sigma_v}.
#' @param sigmae The standard deviation of the observation process \eqn{\sigma_e}.
#' @param x0 The initial state.
#' @param P0 The initial covariance of the state. 
#'
#' @return
#' The function returns a list with the elements:
#' \itemize{
#' \item{xh: The estimate of the filtered state at time \eqn{t=1,...,T}.}
#' \item{ll: The estimate of the log-likelihood.}
#' }
#' @references 
#' Dahlin, J. & Schoen, T. B. "Getting started with particle 
#' Metropolis-Hastings for inference in nonlinear dynamical models." 
#' pre-print, arXiv:1511.01707, 2015.
#' @author 
#' Johan Dahlin <johan.dahlin@liu.se>
#' @note 
#' See Section 4 in the reference for more details.
#' @keywords 
#' ts
#' @export
#' @example ./examples/kf
#' @importFrom stats dnorm

kf <- function(y, phi, sigmav, sigmae, x0, P0) {

  T     <- length(y)
  yhatp <- matrix(x0, nrow=T, ncol=1)
  xhatf <- matrix(x0, nrow=T, ncol=1)
  xhatp <- matrix(x0, nrow=T+1, ncol=1)
  Pp    <- P0
  ll    <- 0
  
  # Set parameters 
  A <- phi
  C <- 1
  Q <- sigmav^2
  R <- sigmae^2
  
  for (tt in 2:T) {
    
    # Compute Kalman Gain
    S <- C * Pp * C + R
    K <- Pp * C / S
    
    # Compute state estimate
    yhatp[tt]   <- C * xhatp[tt]
    xhatf[tt]   <- xhatp[tt] + K * (y[tt] - yhatp[tt])
    xhatp[tt+1] <- A * xhatf[tt]
    
    # Update covariance
    Pf <- Pp - K * S * K
    Pp <- A * Pf * A + Q
    
    # Estimate loglikelihood (not in the last iteration, to be able to compare with faPF)
    if ( tt < T ) { 
      ll = ll + dnorm(y[tt], yhatp[tt], sqrt(S), log=TRUE) 
    }
  }
  
  list(xh = xhatf, ll=ll)
}

##############################################################################
# Bootstrap particle filter (SV model)
##############################################################################

#' Bootstrap particle filter for state estimate in a simple stochastic 
#' volatility model
#' @description 
#' Estimates the filtered state and the log-likelihood for a stochastic 
#' volatility model of the form \eqn{x_t = \mu + \phi ( x_{t-1} - \mu ) + 
#' \sigma_v v_t} and \eqn{y_t = \exp(x_t/2) e_t}, where \eqn{v_t} and \eqn{e_t} 
#' denote independent standard Gaussian random variables, i.e. \eqn{N(0,1)}.
#' @param y Observations from the model for \eqn{t=1,...,T}.
#' @param mu The mean of the log-volatility process \eqn{\mu}.
#' @param phi The persistence of the log-volatility process \eqn{\phi}.
#' @param sigmav The standard deviation of the log-volatility process 
#' \eqn{\sigma_v}.
#' @param nPart The number of particles to use in the filter.
#' @param T The number of observations.
#'
#' @return
#' The function returns a list with the elements:
#' \itemize{
#' \item{xh: The estimate of the filtered state at time \eqn{t=1,...,T}.}
#' \item{ll: The estimate of the log-likelihood.}
#' }
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
#' @example ./examples/sm_sv
#' @importFrom stats dnorm
#' @importFrom stats rnorm

sm_sv <- function(y, mu, phi, sigmav, nPart, T) {

  #===========================================================
  # Initialise variables
  #===========================================================
  a     <- matrix(0, nrow=nPart, ncol=T+1)
  p     <- matrix(0, nrow=nPart, ncol=T+1)
  w     <- matrix(1/nPart, nrow=nPart, ncol=T+1)
  v     <- matrix(1 , nrow=nPart, ncol=T+1)
  ll    <- 0
  
  # Generate initial state
  p[,1] <- rnorm(nPart, mu, sigmav / sqrt( 1 - phi * phi ))
  a[,1] <- 1:nPart
  
  #===========================================================
  # Run main loop
  #===========================================================
  for (tt in 2:(T+1)) {
    
    #=========================================================
    # Resample ( multinomial )
    #=========================================================
    idx <- sample(1:nPart, nPart, replace=TRUE, prob = w[,tt-1])
    
    # Resample the ancestory linage
    a[,1:tt-1]  <- a[idx,1:tt-1]
    
    # Add the most recent ancestors
    a[,tt]      <- idx
    
    #=========================================================
    # Propagate
    #=========================================================
    p[,tt] <- mu + phi * (p[idx,tt-1] - mu) + rnorm(nPart, 0, sigmav) 
    
    #=========================================================
    # Compute weights
    #=========================================================
    v[,tt] <- dnorm(y[tt-1], 0, exp( p[,tt] / 2 ), log=TRUE)
    
    # Rescale log-weights and recover weight
    vmax   <- max(v[,tt])
    v[,tt] <- exp(v[,tt] - vmax)
    
    # Normalize the weights
    w[,tt] <- v[,tt] / sum(v[,tt])
    
    # Estimate the log-likelihood
    ll     <- ll + vmax + log(sum(v[,tt])) - log(nPart)
    
  }
  #===========================================================
  # Return state estimate and log-likelihood estimate
  #===========================================================

  # Sample the state estimate using the weights at tt=T
  nIdx  <- sample(1:nPart, 1, prob=w[,T])
  xhatf <- p[ cbind(a[nIdx,], 1:(T+1)) ]  
  
  list(xh = xhatf, ll=ll)
}

##############################################################################
# End of file
##############################################################################
