#' Steps of the EM algorithm for a Gaussian Mixture model
#'
#' Functions to perform the expectation and maximization steps of the EM
#' algorithm for an multivariate Gaussian mixture model.
#'
#' @aliases EStep MStep
#' @param x A matrix of observations where rows corresponds to features and
#'   columns to experiments.
#' @param theta A list of parameters formatted as described in
#'   \code{\link{rtheta}}.
#' @param kappa A matrix where the (i,j)'th entry is the probability that
#'   \code{x[i,]} belongs to the \code{j}'th component. Usually the returned
#' value of \code{EStep}.
#' @param meta.special.case Logical. If \code{TRUE}, the maximization step is
#'   performed under the special case of Li et. al. (2011). Default values is
#' \code{FALSE}.
#' @return \code{EStep} returns a matrix of probabilities as \code{kappa}
#'   above.
#' @author Anders Ellern Bilgrau <anders.ellern.bilgrau@@gmail.com>
#' @seealso \code{\link{rtheta}}
#' @references Li, Q., Brown, J. B. J. B., Huang, H., & Bickel, P. J. (2011).
#'   Measuring reproducibility of high-throughput experiments. The Annals of
#'   Applied Statistics, 5(3), 1752-1779. doi:10.1214/11-AOAS466
#' @examples
#' set.seed(1)
#' sim <- GMCM:::SimulateGMMData(n = 100)
#' x <- sim$z
#' true.theta <- sim$theta
#' init.theta <- GMCM:::rtheta()  # Generate starting parameters
#'
#' # Do one EM interation
#' es <- GMCM:::EStep(x, init.theta)
#' new.theta <- GMCM:::MStep(x, es)
#'
#' # Compare current estimate with the true
#' new.theta
#' true.theta
#' @keywords internal
EStep <- function (x, theta) {
  return(EStepRcpp(x, mus = theta$mu, sigmas = theta$sigma, pie = theta$pie))
}

# # Old EStep:
# EStep <- function (x, theta) {
#   g <- function(j) {
#     theta$pie[j]*dmvnormal(x, mu = theta$mu[[j]], sigma = theta$sigma[[j]])
#   }
#   f <- sapply(1:theta$m, g)
#   kappa <- f/rowSums(f)
#   kappa[is.nan(kappa[,1]) | is.nan(kappa[,2]), ] <- 0  # What if kappa has more than 3 cols, fix!
#   return(kappa)
# }
