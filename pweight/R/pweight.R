#' pweight: P-value weighting in multiple testing.
#'
#' The pweight package provides functions to compute various p-value
#' weighting schemes, as Spjotvoll, exponential and Bayes weights. These are methods for improving power in multiple testing via the use of prior informaton.
#'
#'
#' Several of these methods were proposed by the authors in the following papers:
#' "Optimal Multiple Testing Under a
#'  Gaussian Prior on the Effect Sizes", by Dobriban, Fortney, Kim and Owen,
#'   \url{http://arxiv.org/abs/1504.02935}
#'
#' @section Weighting functions:
#' \code{\link{bayes_weights}} computes Bayes p-value weights
#'
#' \code{\link{spjotvoll_weights}} computes the Spjotvoll p-value weights
#'
#' \code{\link{exp_weights}} computes the exponential weights
#'
#' @docType package
#' @name pweight
#' @examples
#' J <- 100
#' mu <- rnorm(J)
#' sigma <- 1 * rep(1, J)
#' q <- 0.05 / J
#' res <- bayes_weights(mu, sigma, q)
NULL
