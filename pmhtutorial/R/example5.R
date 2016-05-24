##############################################################################
#
# Example of particle Metropolis-Hastings 
# in a stochastic volatility model
#
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

#' Parameter estimation in a simple stochastic volatility model
#' 
#' @description
#' Minimal working example of parameter estimation in a stochastic volatility 
#' model using the particle Metropolis-Hastings algorithm with a bootstrap 
#' particle filter for providing an unbiased estimator of the likelihood. The 
#' code estimates the parameter posterior for three parameters using 
#' real-world data. 
#' @details 
#' The Particle Metropolis-Hastings (PMH) algorithm makes use of a Gaussian 
#' random walk as the proposal for the parameters. The data are scaled 
#' log-returns from the OMXS30 index during the period from January 2, 2012 
#' to January 2, 2014.
#' 
#' This version of the code makes use of a proposal that is tuned using a 
#' pilot run. Furthermore the model is reparameterised to enjoy better mixing 
#' properties by making the parameters unrestricted to a certain part of the 
#' real-line.
#' @param nIter The number of iterations in the PMH algorithm. 100 iterations 
#' takes about a minute on a laptop to execute. 7500 iterations are used 
#' in the reference below. The length of the burn-in is calculated as one 
#' fifth of nIter.
#' @return 
#' The function returns the estimated marginal parameter posteriors for each 
#' parameter, the trace of the Markov chain and the resulting autocorrelation 
#' function. The data is also presented with an estimate of the 
#' log-volatility.
#' 
#' The function returns a list with the elements:
#' \itemize{
#' \item{thhat: The estimate of the mean of the parameter posterior.}
#' \item{xhat: The estimate of the mean of the log-volatility posterior.}
#' \item{thhatSD: The estimate of the standard deviation of the parameter 
#' posterior.}
#' \item{xhatSD: The estimate of the standard deviation of the log-volatility 
#' posterior.}
#' \item{iact: The estimate of the integrated autocorrelation time for each 
#' parameter.}
#' \item{estCov: The estimate of the covariance of the parameter posterior.}
#' }
#' @references 
#' Dahlin, J. & Schoen, T. B. "Getting started with particle 
#' Metropolis-Hastings for inference in nonlinear dynamical models." 
#' pre-print, arXiv:1511.01707, 2015.
#' @author 
#' Johan Dahlin <johan.dahlin@liu.se>
#' @note 
#' See Section 6.3.2 in the reference for more details.
#' @example ./examples/example5
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

example5_sv <- function( nIter=7500 ) {

  # Set the random seed to replicate results in tutorial
  set.seed(10)
  
  ##############################################################################
  # Define the model
  ##############################################################################
  
  # Here, we use the following model
  #
  # x[tt+1] = phi  * x[tt] + sigma   * v[tt]
  # y[tt]   = beta * exp( xt[tt]/2 ) * e[tt]
  #
  # where v[tt] ~ N(0,1) and e[tt] ~ N(0,1)
  
  # Set the number of time steps to simulate
  T      <- 500;
  
  ##############################################################################
  # Load data
  ##############################################################################
  d <- Quandl::Quandl("NASDAQOMX/OMXS30", start_date="2012-01-02", end_date="2014-01-02", type="zoo")
  y <- as.numeric(100 * diff(log(d$"Index Value")))
  
  ##############################################################################
  # Parameter estimation using PMH
  ##############################################################################
  
  # The inital guess of the parameter
  initPar  <- c( 0, 0.9, 0.2)
  
  # No. particles in the particle filter ( choose nPart ~ T )
  nPart    <- 500
  
  # The length of the burn-in and the no. iterations of PMH ( nBurnIn < nIter )
  nBurnIn  <- floor(nIter / 5)
  
  # The standard deviation in the random walk proposal
  stepSize <- matrix( c( 0.104402823,  0.008707826, -0.009090245, 
                         0.008707826,  0.071040325, -0.034589832,
                        -0.009090245, -0.034589832,  0.086980437), ncol=3, nrow=3)
  stepSize <- 0.8 * stepSize
  
  # Run the PMH algorithm
  res <- pmh_sv_reparameterised(y, initPar, nPart, T, nIter, stepSize)
  
  ##############################################################################
  # Plot the results
  ##############################################################################
  
  # Extract the states after burn-in
  resTh <- res$thhat[nBurnIn:nIter,]
  resXh <- res$xhat[nBurnIn:nIter,]
  
  # Estimate the KDE of the marginal posteriors
  kde1  <- density(resTh[,1], kernel="e", from=-1, to=1)
  kde2  <- density(resTh[,2], kernel="e", to=1)
  kde3  <- density(resTh[,3], kernel="e", from=0)
  
  # Estimate the posterior mean and the corresponding standard deviation
  thhat   <- colMeans(resTh)
  thhatSD <- apply(resTh, 2, sd)
  
  # Estimate the log-volatility and the corresponding standad deviation
  xhat    <- colMeans(resXh)
  xhatSD  <- apply(resXh, 2, sd)
  
  # Plot the parameter posterior estimate, solid black line indicate posterior mean
  # Plot the trace of the Markov chain after burn-in, solid black line indicate posterior mean
  layout(matrix(c( 1, 1, 1, 2, 2, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11), 5, 3, byrow=TRUE))
  par(mar=c(4, 5, 0, 0))
  
  # Grid for plotting the data and log-volatility
  gridy <- seq(1, length(y), 1)
  
  plot(y, col="#1B9E77", lwd=1, type="l", xlab="time", ylab="log-returns", ylim=c(-5,5), bty="n")
  plot(xhat[-1], col="#D95F02", lwd=1.5, type="l", xlab="time", ylab="log-volatility estimate", ylim=c(-2,2), bty="n")
  
  nPlot <- nIter - nBurnIn
  grid  <- seq(nBurnIn, nBurnIn + nPlot - 1, 1)
  
  # Mu
  hist( resTh[,1], breaks=floor(sqrt(nIter - nBurnIn)), col=rgb(t(col2rgb("#7570B3")) / 256, alpha=0.25), border=NA,
        xlab=expression(mu), ylab="posterior estimate", main="", xlim=c(-1,1), freq=FALSE)
  lines(kde1, lwd=2, col="#7570B3" ); 
  abline(v=thhat[1], lwd=1,lty="dotted" );
  
  plot(resTh[,1], col='#7570B3', type="l", xlab="iteration (after burn-in)", ylab=expression(mu), ylim=c(-1,1), bty="n")
  abline(h=thhat[1], lwd=1, lty="dotted")
  
  muACF <- acf(resTh[,1], plot=FALSE, lag.max=100)
  plot(muACF$lag, muACF$acf, col = '#7570B3', type="l", xlab="iteration", ylab=expression("ACF of "* mu), lwd=2, ylim=c(-0.5,1), bty="n")
  abline(h=1.96/sqrt(nIter-nBurnIn), lty="dotted")
  abline(h=-1.96/sqrt(nIter-nBurnIn), lty="dotted")
  
  # Phi
  hist( resTh[,2], breaks=floor(sqrt(nIter - nBurnIn)), col=rgb(t(col2rgb("#E7298A")) / 256, alpha=0.25), border=NA,
        xlab=expression(phi), ylab="posterior estimate", main="", xlim=c(0.88,1.0), freq=FALSE)
  lines(kde2, lwd=2, col="#E7298A" ); 
  abline(v=thhat[2], lwd=1,lty="dotted" );
  
  plot(resTh[,2], col='#E7298A', type="l", xlab="iteration (after burn-in)", ylab=expression(phi), ylim=c(0.88,1.0), bty="n")
  abline(h=thhat[2], lwd=1, lty="dotted")
  
  phiACF <- acf(resTh[,2], plot=FALSE, lag.max=100)
  plot(phiACF$lag, phiACF$acf, col = '#E7298A', type="l", xlab="iteration", ylab=expression("ACF of "* phi), lwd=2, ylim=c(-0.5,1), bty="n")
  abline(h=1.96/sqrt(nIter-nBurnIn), lty="dotted")
  abline(h=-1.96/sqrt(nIter-nBurnIn), lty="dotted")
  
  # Sigma[v]
  hist( resTh[,3], breaks=floor(sqrt(nIter - nBurnIn)), col=rgb(t(col2rgb("#66A61E")) / 256, alpha=0.25), border=NA,
        xlab=expression(sigma[v]), ylab="posterior estimate", main="", xlim=c(0.0,0.4), freq=FALSE)
  lines(kde2, lwd=2, col="#66A61E" ); 
  abline(v=thhat[3], lwd=1,lty="dotted" );
  
  plot(resTh[,3], col='#66A61E', type="l", xlab="iteration (after burn-in)", ylab=expression(sigma[v]), ylim=c(0.0,0.4), bty="n")
  abline(h=thhat[3], lwd=1, lty="dotted")
  
  sigmavACF <- acf(resTh[,3], plot=FALSE, lag.max=100)
  plot(sigmavACF$lag, sigmavACF$acf, col = '#66A61E', type="l", xlab="iteration", ylab=expression("ACF of "* sigma[v]), lwd=2, ylim=c(-0.5,1), bty="n")
  abline(h=1.96/sqrt(nIter-nBurnIn), lty="dotted")
  abline(h=-1.96/sqrt(nIter-nBurnIn), lty="dotted")
  
  # Compute an estimate of the IACT using the first 100 ACF coefficients
  iact   <- 1 + 2 * c(sum(muACF$acf), sum(phiACF$acf), sum(sigmavACF$acf))
  
  # Estimate the covariance of the posterior to tune the proposal
  estCov <- var( resTh )

  # Compile output
  list(thhat=thhat, xhat=xhat, thhatSD=thhatSD, xhatSD=xhatSD, iact=iact, estCov=estCov)
}

##############################################################################
# End of file
##############################################################################
