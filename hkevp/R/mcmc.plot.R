#' @title
#' Markov chains plotting
#' 
#' @author
#' Quentin Sebille
#' 
#' 
#' 
#' @description 
#' Plot of the resulting Markov chains obtained by the MCMC procedure \code{hkevp.fit}. May be used to assess graphically convergence of the chains.
#'
#'
#'
#'
#' 
#' @param fit
#' A named list.
#' Output from the \code{hkevp.fit} procedure.
#' 
#' @param mfrow
#' A vector of two numerical values.
#' Parameter of the window plotting called by the \code{plot(...)} function. Optional.
#' 
#' @param plot.spatial
#' Logical.
#' Should the Markov chains of the sills and ranges hyperparameters should be plotted (FALSE by default)?
#'
#' 
#' 
#' 
#' @seealso \link{mcmc.accept}
#' 
#'
#' @examples
#' # Simulation of HKEVP:
#' sites <- as.matrix(expand.grid(1:3,1:3))
#' knots <- sites
#' mu <- sites[,1]*10
#' sigma <- 3
#' xi <- .2
#' alpha <- .4
#' tau <- 1
#' ysim <- hkevp.rand(10, sites, knots, mu, sigma, xi, alpha, tau)
#' 
#' # HKEVP fit:
#' fit <- hkevp.fit(ysim, sites, knots, niter = 100, nburn = 50, quiet = FALSE)
#' 
#' # Markov chains plot:
#' # mcmc.plot(fit)
#' 
#' 
#' 
mcmc.plot <- function(fit, mfrow, plot.spatial) {
  # Default value
  if (missing(plot.spatial)) plot.spatial <- FALSE
  if (missing(mfrow))
    if (plot.spatial) {mfrow <- c(3,3)} else {mfrow <- c(2,3)}
  
  # Plotting parameters
  nsites <- dim(fit$GEV)[1]
  COLORS <- rgb(0:(nsites - 1)/nsites,0:(nsites - 1)/nsites,0:(nsites - 1)/nsites)
  par(mfrow = mfrow)
  
  
  # 1/ GEV parameters
  if (fit$spatial$vary[1])
    matplot(t(fit$GEV[,1,]), type = 'l', lty = 1, col = COLORS, xlab = 'Iterations', ylab = '', axes = FALSE)
  else 
    plot(fit$GEV[1,1,], type = 'l', lty = 1, col = COLORS, xlab = 'Iterations', ylab = '', axes = FALSE)
  box()
  axis(1)
  axis(2, las = 2)
  title(expression(mu))
  
  if (fit$spatial$vary[2])
    matplot(t(fit$GEV[,2,]), type = 'l', lty = 1, col = COLORS, xlab = 'Iterations', ylab = '', axes = FALSE)
  else 
    plot(fit$GEV[1,2,], type = 'l', lty = 1, col = COLORS, xlab = 'Iterations', ylab = '', axes = FALSE)
  box()
  axis(1)
  axis(2, las = 2)
  title(expression(sigma))
  
  if (fit$spatial$vary[3])
    matplot(t(fit$GEV[,3,]), type = 'l', lty = 1, col = COLORS, xlab = 'Iterations', ylab = '', axes = FALSE)
  else 
    plot(fit$GEV[1,3,], type = 'l', lty = 1, col = COLORS, xlab = 'Iterations', ylab = '', axes = FALSE)
  box()
  axis(1)
  axis(2, las = 2)
  title(expression(xi))
  
  
  # 2/ Dependence parameters
  plot(fit$alpha, type = 'l', xlab = 'Iterations', ylab = '', axes = FALSE)
  box()
  axis(1)
  axis(2, las = 2)
  title(expression(alpha))
  
  plot(fit$tau, type = 'l', xlab = 'Iterations', ylab = '', axes = FALSE)
  box()
  axis(1)
  axis(2, las = 2)
  title(expression(tau))
  
  
  # 3/ Log-likelihood
  plot(fit$llik, type = 'l', xlab = 'Iterations', ylab = '', axes = FALSE)
  box()
  axis(1)
  axis(2, las = 2)
  title(expression("log-likelihood"))
  
  
  # 4/ GEV Spatial parameters (if plot.spatial is TRUE, optional)
  if (plot.spatial) {
    # Sills
    matplot(fit$spatial$sills, type = 'l', lty = 1, col = fit$spatial$vary, xlab = 'Iterations', ylab = '', axes = FALSE)
    box()
    axis(1)
    axis(2, las = 2)
    title(expression("Sills "*delta))
    
    matplot(fit$spatial$ranges, type = 'l', lty = 1, col = fit$spatial$vary, xlab = 'Iterations', ylab = '', axes = FALSE)
    box()
    axis(1)
    axis(2, las = 2)
    title(expression("Ranges "*lambda))
    
  }
  
}
