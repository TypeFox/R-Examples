#' EM algorithm for Gaussian mixture models
#'
#' The regular expectation-maximization algorithm for general multivariate
#' Gaussian mixture models.
#'
#' Though not as versatile, the algorithm can be a faster alternative to
#' \code{Mclust} in the \code{mclust}-package.
#'
#' @param x A matrix of observations where rows correspond to features and
#'   columns to experiments.
#' @param theta A list of parameters as described in \code{\link{rtheta}}.
#' @param eps The maximal required difference in successive likelihoods to
#'   establish convergence.
#' @param max.ite The maximum number of iterations.
#' @param trace.theta Logical. If \code{TRUE}, all estimates are stored and
#'   returned. Default is \code{FALSE}.
#' @param verbose Set to \code{TRUE} for verbose output. Default is
#' \code{FALSE}.
#' @return
#'   A list of length 3 with elements:
#'   \item{theta }{A list of the estimated parameters as described in
#'                 \code{\link{rtheta}}.}
#'   \item{loglik.tr}{A numeric vector of the log-likelihood trace.}
#'   \item{kappa }{A matrix where \code{kappa[i,j]} is the probability that
#'                 \code{x[i, ]} is realized from the \code{j}'th component.}
#' @author Anders Ellern Bilgrau <anders.ellern.bilgrau@@gmail.com>
#' @seealso \code{\link{rtheta}}, \code{\link{PseudoEMAlgorithm}}
#' @examples
#'
#' set.seed(10)
#' data <- SimulateGMCMData(n = 1000, d = 2, m = 3)
#' start.theta <- rtheta(d = 2, m = 3)
#' res <- GMCM:::EMAlgorithm(data$z, theta = start.theta)
#'
#' par(mfrow = c(1,2))
#' plot(data$z, cex = 0.5, pch = 16, main = "Simulated data",
#'      col = rainbow(3)[data$K])
#' plot(data$z, cex = 0.5, pch = 16, main = "GMM clustering",
#'      col = rainbow(3)[apply(res$kappa,1,which.max)])
#'
EMAlgorithm <- function (x, theta, eps = 1e-6, max.ite = 1e5,
                         trace.theta = FALSE, verbose = FALSE) {
  loglik.tr <- c(dgmm.loglik(theta, x))
  theta.tr  <- vector("list", 1); theta.tr[[1]] <- theta
  for (k in 2:max.ite) {
    kappa <- EStep(x = x, theta = theta)
    theta <- MStep(x = x, kappa = kappa)
    loglik.tr[k] <- dgmm.loglik(theta, x)
    theta.tr[[k]] <- theta
    delta <- loglik.tr[k] - loglik.tr[k-1]
    if (verbose) {
      cat("iteration", k, "\tDelta loglik =", delta, "\n"); flush.console()
    }
    if (delta < 0)
      stop("Delta likelihood was negative. Something went wrong!")
    if (delta < eps)
      break
    if (k == max.ite)
      warning(paste("Max (", max.ite, ") iterations reached", sep = ""))
  }
  res <- list(theta = theta,
              loglik.tr = loglik.tr,
              kappa = kappa,
              theta.tr = theta.tr)
  if (!trace.theta)
    res <- res[-4]
  return(res)
}
