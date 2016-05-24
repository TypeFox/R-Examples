#' @title
#' Acceptance ratio of the HKEVP fit
#' 
#' @author
#' Quentin Sebille
#' 
#' 
#' 
#' @description
#' Computation and plot of the acceptance rates for each parameter of the HKEVP, from the output of the MCMC procedure.
#'
#'
#'
#'
#'
#' @param fit
#' A named list.
#' Output from the \code{hkevp.fit} procedure.
#' 
#' @param plotted
#' Logical.
#' Whether or not the acceptance rates should be plotted.
#'
#'
#'
#' @return 
#' A named list with four elements:
#' \itemize{
#' \item \code{GEV}: A matrix of real values. Acceptance rates for each GEV parameter (columns) and each site position (rows).
#' \item \code{alpha}: A numerical value. Acceptance ratio for the dependence parameter \eqn{\alpha}.
#' \item \code{tau}: A numerical value. Acceptance ratio for the bandwidth parameter \eqn{\tau}.
#' \item \code{A}: A matrix of real values. Acceptance rates of the random effect for each knot position (rows) and each block (columns) where maxima are observed.
#' }
#' 
#' 
#' @seealso \link{mcmc.plot}
#' 
#' 
#' @export
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
#' # Acceptance rates:
#' # mcmc.accept(fit, TRUE)
#' 
#' 
#' 
#' 
mcmc.accept <- function(fit, plotted) {
  # Default value
  if (missing(plotted)) plotted <- TRUE
  
  # Computing the acceptance ratio
  compute.ratio <- function(x) sum(diff(x) != 0) / (length(x) - 1)
  GEV.acc <- apply(fit$GEV, 1:2, compute.ratio)
  alpha.acc <- compute.ratio(fit$alpha)
  tau.acc <- compute.ratio(fit$tau)
  A.acc <- apply(X = fit$A, MARGIN = 1:2, compute.ratio)
  
  # Plotting the acceptance ratio
  if (plotted) {
    par(mfrow = c(1,1))
    suppressWarnings(boxplot(cbind(GEV.acc, alpha.acc, tau.acc, as.vector(A.acc)), ylim = 0:1, axes = FALSE, main = expression("Acceptance rates")))
    box(); axis(2, las = 1); axis(1, at = 1:6, labels = c(expression(mu), expression(sigma), expression(xi), expression(alpha), expression(tau), expression(A)))
  }
  
  result <- list(GEV = GEV.acc, alpha = alpha.acc, tau = tau.acc, A = A.acc)
  
  return(result)  
  
}





