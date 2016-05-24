#' iLaplace: A package for Approximating Multidimensional Integrals of Unimodal Functions
#'
#' This package implements the improved Laplace approximation of Ruli et al. (2015) for multivariate
#' integrals of user-written unimodal functions. This method essentially approximates the target
#' integral by the ratio of the integrand and its normalised version, both evaluated at the modal value.
#' The normalised integrand is obtained through a sequential application of the Laplace approximation
#' for marginal densities. Like the standard Laplace approximation, the improved Laplace approximation
#' is a deterministic method which approximates intractable multidimensional integrals by (essentially)
#' numerical optimisations. However, whith respect to the Laplace approximation, the improved Laplace
#' involves scalar numerical integrations. Nevertheless, the improved Laplace approximation tends to be
#' extremely accurate, especially in those cases in which the integrand is skewed or has fat "tails".
#'
#' @docType package
#' @name iLaplace-package
#' @author Erlis Ruli \email{erlisr@@yahoo.it}
#' @seealso \code{\link[iLaplace]{iLap}} and \code{\link[iLaplace]{iLap2d}} for more details and examples.
#' @references
#' Ruli E., Sartori N. & Ventura L. (2015)
#' Improved Laplace approximation for marignal likelihoods. \url{http://arxiv.org/abs/1502.06440}
NULL
