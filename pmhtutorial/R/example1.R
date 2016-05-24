##############################################################################
#
# Example of fully-adapted particle filtering 
# in a linear Gaussian state space model
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

#' State estimation in a linear Gaussian state space model
#' 
#' @description
#' Minimal working example of state estimation in a linear Gaussian state 
#' space model using Kalman filtering and a fully-adapted particle filter. 
#' The code estimates the bias and mean squared error (compared with the 
#' Kalman estimate) while varying the number of particles in the particle 
#' filter.
#' @details
#' The Kalman filter is a standard implementation without an input. The 
#' particle filter is fully adapted (i.e. takes the current observation into 
#' account when proposing new particles and computing the weights).
#' @return 
#' Returns a plot with the generated observations y and the difference in the
#' state estimates obtained by the Kalman filter (the optimal solution) and 
#' the particle filter (with 20 particles). Furthermore, the function returns 
#' plots of the estimated bias and mean squared error of the state estimate 
#' obtained using the particle filter (while varying the number of particles) 
#' and the Kalman estimates.
#' 
#' The function returns a list with the elements:
#' \itemize{
#' \item{y: The observations generated from the model.}
#' \item{x: The states generated from the model.}
#' \item{xhatKF: The estimate of the state from the Kalman filter.}
#' \item{xhatPF: The estimate of the state from the particle filter with 
#' 20 particles.}
#' }
#' @references 
#' Dahlin, J. & Schoen, T. B. "Getting started with particle 
#' Metropolis-Hastings for inference in nonlinear dynamical models." 
#' pre-print, arXiv:1511.01707, 2015.
#' @author 
#' Johan Dahlin <johan.dahlin@liu.se>
#' @note 
#' See Section 4.2 in the reference for more details.
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

example1_lgss <- function() {

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
  phi     = 0.75
  sigmav  = 1.00
  sigmae  = 0.10
  
  # Set the number of time steps to simulate
  T       = 250
  
  # Set the initial state
  x0      = 0
  
  ##############################################################################
  # Generate data
  ##############################################################################
  
  data <- generateData( phi, sigmav, sigmae, T, x0 )
  x    <- data$x
  y    <- data$y
  
  # Plot the latent state and observations
  layout( matrix( c(1,1,2,2,3,4), 3, 2, byrow = TRUE )  )
  par   ( mar = c(4,5,0,0) )
  
  grid <- seq( 0, T )
  
  plot( grid, y, col="#1B9E77", type="l", xlab="time", ylab="observation", ylim=c(-6,6), bty="n")
  
  ###################################################################################
  # State estimation using the particle filter and Kalman filter
  ###################################################################################
  
  # Using N = 20 particles and plot the estimate of the latent state
  N      <- 20
  outPF  <- sm( y, phi, sigmav, sigmae, N, T, x0 )
  outKF  <- kf( y, phi, sigmav, sigmae, x0, 0.01 )
  diff   <- outPF$xh - outKF$xh[-(T+1)]
  
  grid   <- seq( 0, T-1 )
  plot( grid, diff, col="#7570B3", type="l", xlab="time", ylab="difference in state estimate", ylim=c(-0.1,0.1), bty="n")
  
  # Compute bias and MSE
  logBiasMSE = matrix( 0, nrow = 7, ncol = 2 )
  gridN      = c( 10, 20, 50, 100, 200, 500, 1000)
  
  for ( ii in 1:length(gridN) ) {
    smEstimate          <- sm(y,phi,sigmav,sigmae,gridN[ii],T,x0)$xh
    kmEstimate          <- outKF$xh[-(T+1)]
    
    logBiasMSE[ ii, 1 ] <- log( mean( abs( smEstimate - kmEstimate )   ) )  
    logBiasMSE[ ii, 2 ] <- log( mean(    ( smEstimate - kmEstimate )^2 ) )
  }
  
  # Plot the bias and MSE for comparison
  plot(   gridN, logBiasMSE[ , 1 ], col="#E7298A", type="l", xlab="no. particles (N)", ylab="log-bias", ylim=c(-7,-3), bty="n")
  points( gridN, logBiasMSE[ , 1 ], col="#E7298A", pch=19 )
  
  plot(   gridN, logBiasMSE[ , 2 ], col="#66A61E", lwd=1.5, type="l", xlab="no. particles (N)", ylab="log-MSE", ylim=c(-12,-6), bty="n")
  points( gridN, logBiasMSE[ , 2 ], col="#66A61E", pch=19 )

  list(y=y, x=x, xhatKF=outKF$xh[-(T+1)], xhatPF=outPF$xh)  
}

###################################################################################
# End of file
###################################################################################
