##############################################################################
#
# Example of particle Metropolis-Hastings 
# in a linear Gaussian state space model
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

#' Parameter estimation in a linear Gaussian state space model
#' 
#' @description
#' Minimal working example of parameter estimation in a linear Gaussian state 
#' space model using the particle Metropolis-Hastings algorithm with a 
#' fully-adapted particle filter for providing an unbiased estimator of the 
#' likelihood. The code estimates the parameter posterior for one parameter 
#' using simulated data. 
#' @details 
#' The Particle Metropolis-Hastings (PMH) algorithm makes use of a Gaussian 
#' random walk as the proposal for the parameter. The PMH algorithm is run 
#' using different step lengths in the proposal. This is done to illustrate 
#' the difficulty when tuning the proposal and the impact of a too 
#' small/large step length.
#' @param nIter The number of iterations in the PMH algorithm. 100 iterations 
#' takes about ten seconds on a laptop to execute. 5000 iterations are used 
#' in the reference below. The length of the burn-in is calculated as one 
#' fifth of nIter.
#' @return 
#' Returns the estimate of the posterior mean.
#' @references 
#' Dahlin, J. & Schoen, T. B. "Getting started with particle 
#' Metropolis-Hastings for inference in nonlinear dynamical models." 
#' pre-print, arXiv:1511.01707, 2015.
#' @author 
#' Johan Dahlin <johan.dahlin@liu.se>
#' @note 
#' See Section 4.4 in the reference for more details.
#' @example ./examples/example2
#' @keywords 
#' misc
#' @export
#' @importFrom grDevices col2rgb
#' @importFrom grDevices rgb
#' @importFrom graphics abline
#' @importFrom graphics hist
#' @importFrom graphics layout
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom graphics points
#' @importFrom stats acf
#' @importFrom stats density
#' @importFrom stats sd
#' @importFrom stats var

example2_lgss <- function( nIter=5000 ) {
  
  # Set the random seed to replicate results in tutorial
  set.seed( 10 )
  
  ##############################################################################
  # Define the model
  ##############################################################################
  
  # Here, we use the following model
  #
  # x[tt+1] = phi   * x[tt] + sigmav * v[tt]
  # y[tt]   = x[tt]         + sigmae * e[tt]
  #
  # where v[tt] ~ N(0,1) and e[tt] ~ N(0,1)
  
  # Set the parameters of the model
  phi     <- 0.75
  sigmav  <- 1.00
  sigmae  <- 0.10
  
  # Set the number of time steps to simulate
  T      <- 250
  
  # Set the initial state
  x0     <- 0
  
  
  ##############################################################################
  # Generate data
  ##############################################################################
  
  data <- generateData( phi, sigmav, sigmae, T, x0 )
  
  ##############################################################################
  # Parameter estimation using PMH
  ##############################################################################
  
  # The inital guess of the parameter
  initPar  <- 0.50
  
  # No. particles in the particle filter
  nPart    <- 100
  
  # The length of the burn-in and the no. iterations of PMH ( nBurnIn < nIter )
  nBurnIn  <- floor(nIter / 5)
  
  # Run the PMH algorithm
  res1 <- pmh(data$y, initPar, sigmav, sigmae, nPart, T, x0, nIter, stepSize = 0.01)
  res2 <- pmh(data$y, initPar, sigmav, sigmae, nPart, T, x0, nIter, stepSize = 0.10)
  res3 <- pmh(data$y, initPar, sigmav, sigmae, nPart, T, x0, nIter, stepSize = 0.50)
  
  ##############################################################################
  # Plot the results
  ##############################################################################
  
  layout( matrix(1:9, 3, 3, byrow = TRUE) ) 
  par   ( mar = c(4,5,0,0) )
  
  # Plot the parameter posterior estimate
  hist(res1[nBurnIn:nIter], breaks=floor(sqrt(nIter-nBurnIn)), col=rgb(t(col2rgb("#7570B3"))/256, alpha=0.25), border=NA, 
       xlab=expression(phi), ylab="posterior estimate", main="", xlim=c(0.5,0.8), ylim=c(0,12), freq=FALSE)
  abline(v=mean(res1[nBurnIn:nIter]), lwd=1, lty="dotted")
  
  hist(res2[nBurnIn:nIter], breaks=floor(sqrt(nIter-nBurnIn)), col=rgb(t(col2rgb("#E7298A"))/256, alpha=0.25), border=NA, 
       xlab=expression(phi), ylab="posterior estimate", main="", xlim=c(0.5,0.8), ylim=c(0,12), freq=FALSE)
  abline(v=mean(res2[nBurnIn:nIter]), lwd=1, lty="dotted" )
  
  hist(res3[nBurnIn:nIter], breaks=floor(sqrt(nIter-nBurnIn)), col=rgb(t(col2rgb("#66A61E"))/256, alpha=0.25), border=NA, 
       xlab=expression(phi), ylab="posterior estimate", main="", xlim=c(0.5,0.8), ylim=c(0,12), freq=FALSE)
  abline(v=mean(res3[nBurnIn:nIter]), lwd=1, lty="dotted" )
  
  # Plot the trace of the Markov chain after the burn-in
  grid <- seq(nBurnIn, nIter, 1)
  
  plot(grid, res1[grid], col='#7570B3', type="l", xlab="iteration (after burn-in)", ylab=expression(phi), 
       ylim=c(0.5,0.8), bty="n")
  abline(h=mean(res1[grid]), lwd=1, lty="dotted")
  
  plot(grid, res2[grid], col='#E7298A', type="l", xlab="iteration (after burn-in)", ylab=expression(phi), 
       ylim=c(0.5,0.8), bty="n")
  abline(h=mean(res2[grid]), lwd=1, lty="dotted")
  
  plot(grid, res3[grid], col='#66A61E', type="l", xlab="iteration (after burn-in)", ylab=expression(phi), 
       ylim=c(0.5,0.8), bty="n")
  abline(h=mean(res3[grid]), lwd=1, lty="dotted")
  
  # Plot the ACF of the Markov chain
  grid <- seq(nBurnIn, nIter, 1)
  
  res1ACF <- acf(res1[grid], plot=FALSE, lag.max=60)
  plot(res1ACF$lag, res1ACF$acf, col='#7570B3', type="l", xlab="iteration", ylab="ACF", ylim=c(-0.2,1), bty="n")
  abline(h=1.96/sqrt(length(grid)), lty="dotted")
  abline(h=-1.96/sqrt(length(grid)), lty="dotted")
  
  res2ACF <- acf(res2[grid], plot=FALSE, lag.max=60)
  plot(res2ACF$lag, res2ACF$acf, col='#E7298A', type="l", xlab="iteration", ylab="ACF", ylim=c(-0.2,1), bty="n")
  abline(h=1.96/sqrt(length(grid)), lty="dotted")
  abline(h=-1.96/sqrt(length(grid)), lty="dotted")
  
  res3ACF <- acf(res3[grid], plot=FALSE, lag.max=60)
  plot(res3ACF$lag, res3ACF$acf, col='#66A61E', type="l", xlab="iteration", ylab="ACF", ylim=c(-0.2,1), bty="n")
  abline(h=1.96/sqrt(length(grid)), lty="dotted")
  abline(h=-1.96/sqrt(length(grid)), lty="dotted")
  
  # Estimate the parameter posterior mean
  c(mean(res1[grid]), mean(res2[grid]), mean(res3[grid]))
}

##############################################################################
# End of file
##############################################################################
