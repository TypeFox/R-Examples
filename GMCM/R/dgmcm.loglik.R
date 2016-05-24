#' Probability, density, and likelihood functions of the Gaussian mixture
#' (copula) model
#'
#' Marginal and simultaneous cumulative distribution, log probability density,
#' and log-likelihood functions of the Gaussian mixture model (GMM) and
#' Gaussian mixture copula model (GMCM) and the relevant inverse marginal
#' quantile functions.
#'
#' \code{qgmm.marginal} distributes approximately \code{res} points around the
#' cluster centers according to the mixture proportions in \code{theta$pie} and
#' evaluates \code{pgmm.marginal} on these points. An approximate inverse of
#' \code{pgmm.marginal} function is constructed by linear interpolation of the
#' flipped evaluated coordinates.
#'
#' @aliases dgmcm.loglik dgmm.loglik dgmm.loglik.marginal pgmm.marginal
#'   qgmm.marginal
#' @param theta A list parameters as described in \code{\link{rtheta}}.
#' @param z A matrix of realizations from the latent process where each row
#'   corresponds to an observation.
#' @param u A matrix of (estimates of) realizations from the GMCM where each
#'   row corresponds to an observation.
#' @param marginal.loglik Logical. If \code{TRUE}, the marginal log-likelihood
#'   functions for each multivariate observation (i.e. the log densities) are
#'   returned. In other words, if \code{TRUE} the sum of the marginal
#'   likelihoods is not computed.
#' @param x A matrix where each row corresponds to an observation.
#' @param res The resolution at which the inversion of \code{qgmm.marginal} is
#'   done. Default is 1000.
#' @param spread The number of marginal standard deviations from the marginal
#'   means the \code{pgmm.marginal} is to be evaluated on.
#' @param \dots Arguments passed to \code{qgmm.marginal}.
#' @return The returned value depends on the value of \code{marginal.loglik}.
#'   If \code{TRUE}, the non-summed marginal likelihood values are returned. If
#'   \code{FALSE}, the scalar sum log-likelihood is returned.
#'
#'   \code{dgmcm.loglik}: As above, with the GMCM density.
#'
#'   \code{dgmm.loglik}: As above, with the GMM density.
#'
#'   \code{dgmm.loglik.marginal}: As above, where the j'th element is evaluated
#'   in the j'th marginal GMM density.
#'
#'   \code{pgmm.marginal}: A matrix where the (i,j)'th entry is the (i,j)'th
#'   entry of \code{z} evaluated in the jth marginal GMM density.
#'
#'   \code{qgmm.marginal}: A matrix where the (i,j)'th entry is the (i,j)'th
#'   entry of \code{u} evaluated in the inverse jth marginal GMM density.
#' @author Anders Ellern Bilgrau <anders.ellern.bilgrau@@gmail.com>
#' @examples
#' set.seed(1)
#' data <- SimulateGMCMData(n = 10)
#' u <- data$u
#' z <- data$z
#' print(theta <- data$theta)
#'
#' GMCM:::dgmcm.loglik(theta, u, marginal.loglik = FALSE)
#' GMCM:::dgmcm.loglik(theta, u, marginal.loglik = TRUE)
#'
#' GMCM:::dgmm.loglik(theta, z, marginal.loglik = FALSE)
#' GMCM:::dgmm.loglik(theta, z, marginal.loglik = TRUE)
#'
#' GMCM:::dgmm.loglik.marginal(theta, z, marginal.loglik = FALSE)
#' GMCM:::dgmm.loglik.marginal(theta, z, marginal.loglik = TRUE)
#'
#' GMCM:::pgmm.marginal(z, theta)
#' GMCM:::qgmm.marginal(u, theta)
#' @keywords internal
dgmcm.loglik <- function (theta, u, marginal.loglik = FALSE, ...) {
  z      <- qgmm.marginal(rbind(u), theta, ...)
  tmp    <- dgmm_loglik_marginal(mus = theta$mu, sigmas = theta$sigma,
                                 pie = theta$pie, z = z,
                                 marginal_loglik = marginal.loglik)
  loglik <- dgmm_loglik(mus = theta$mu, sigmas = theta$sigma, pie = theta$pie,
                        z = z, marginal_loglik = marginal.loglik) - rowSums(tmp)
  if (!marginal.loglik)
    loglik <- sum(loglik)
  return(loglik)
}

#
# dgmcm.loglik2 <- function (theta, u, marginal.loglik = FALSE, ...) {
#  z      <- qgmm.marginal(rbind(u), theta, ...)
#  tmp    <- dgmm.loglik.marginal(theta, z)
#  loglik <- dgmm.loglik(theta, z, marginal.loglik = TRUE) - rowSums(tmp)
#  if (!marginal.loglik)
#    loglik <- sum(loglik)
#  return(loglik)
# }
#
# dgmcm.loglik3 <- function (theta, u, marginal.loglik = FALSE, ...) {
#   am <- function(u) {
#     if (is.null(dim(u)))
#       dim(u) <- c(1, length(u))
#     return(u)
#   }
#   z      <- am(qgmm.marginal(am(u), theta, ...))
#   tmp    <- am(dgmm.loglik.marginal(theta, z))
#   loglik <- dgmm.loglik(theta, z, marginal.loglik = TRUE) - rowSums(tmp)
#   if (!marginal.loglik)
#     loglik <- sum(loglik)
#   return(loglik)
# }
#

