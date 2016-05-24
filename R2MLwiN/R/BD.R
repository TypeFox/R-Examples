#' Calculates Brooks-Draper diagnostic
#' 
#' An internal function, for use in \code{\link{sixway}}, which calculates the
#' Brooks-Draper diagnostic, based on an unpublished paper by David Draper. It
#' estimates the length of a Markov chain required to produce a mean estimate
#' to k significant figures with a given accuracy (alpha). See Browne (2012)
#' for further details.
#' 
#' @param est Numeric scalar for the mean of the distribution
#' @param var Numeric scalar for the variance of the distribution
#' @param rho The first lag (i.e. after zero) of the auto-correlation function (ACF) diagnostic
#' @param k Integer scalar corresponding to the number of significant figures (defaults to \code{2})
#' @param alpha Numeric scalar indicating the desired accuracy (defaults to \code{0.05})
#' 
#' @return The Brooks-Draper diagnostic statistic is returned.
#' 
#' @references
#' Browne, W.J. (2012) MCMC Estimation in MLwiN, v2.26.
#' Centre for Multilevel Modelling, University of Bristol.
#' 
#' @author Zhang, Z., Charlton, C.M.J., Parker, R.M.A., Leckie, G., and Browne,
#' W.J. (2015) Centre for Multilevel Modelling, University of Bristol.
#'
#' @seealso
#' \code{\link{sixway}}
#'
BD <- function(est, var, rho, k = 2, alpha = 0.05) {
  # Based on an upublished paper by David Draper
  ceiling(4 * qnorm(1 - alpha/2)^2 * (sqrt(var)/(10^(floor(log10(abs(est))) - k + 1)))^2 * (1 + rho)/(1 - rho))
} 
