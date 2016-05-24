#' @title List of common arguments
#' @name common-arguments
#' @description
#' Common arguments used in several functions in this package.
#' 
#' @param series a \eqn{T \times K} array with \code{T} observations from the 
#' \eqn{K}-dimensional time series \eqn{\mathbf{X}_t}. Can be a \code{matrix}, \code{data.frame}, 
#' or a multivariate \code{ts} object.
#' @param U a \eqn{T \times K} array with \code{T} observations from the 
#' \eqn{K}-dimensional \strong{whitened} (\code{\link{whiten}}) 
#' time series \eqn{\mathbf{U}_t}. Can be a \code{matrix}, \code{data.frame}, or a 
#' multivariate \code{ts} object.
#' @param mvspectrum.output an object of class \code{"mvspectrum"} representing
#' the multivariate spectrum of \eqn{\mathbf{X}_t} (not necessarily \code{normalize}d).
#' @param f.U multivariate spectrum of class \code{'mvspectrum'} with 
#' \code{normalize = TRUE}. 
#  Must add up to \eqn{0.5 \times I_K}, where \eqn{I_K} 
#  is the \eqn{K \times K} identity matrix.
#' @param algorithm.control list; control settings for any \emph{iterative} ForeCA 
#' algorithm. See \code{\link{complete_algorithm_control}} for details.
#' @param entropy.control list; control settings for entropy estimation.
#' See \code{\link{complete_entropy_control}} for details.
#' @param spectrum.control list; control settings for spectrum estimation. 
#' See \code{\link{complete_spectrum_control}} for details.
#' @param entropy.method string; method to estimate the entropy from discrete
#' probabilities \eqn{p_i}; here \emph{probabilities} are the spectral density
#' evaluated at the Fourier frequencies, 
#' \eqn{\widehat{p}_i = \widehat{f}(\omega_i)}.
#' @param spectrum.method string; method for spectrum estimation; see \code{method}
#' argument in \code{\link{mvspectrum}}.
#' @param threshold numeric; values of spectral density below \code{threshold} are set to
#' \eqn{0}; default \code{threshold = 0}.
#' @param smoothing logical; if \code{TRUE} the spectrum will be
#' smoothed with a nonparametric estimate using \code{\link[mgcv]{gam}} 
#' and an exponential family (with \code{link = log}). Only works
#' for univariate spectrum. The smoothing
#' parameter is chosen automatically using generalized cross-validation 
#' (see \code{\link[mgcv]{gam}} for details). Default: \code{FALSE}.
#' @param base logarithm base; entropy is measured in ``nats'' for 
#' \code{base = exp(1)}; in ``bits'' if \code{base = 2} (default).
NULL