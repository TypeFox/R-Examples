#' Simulation from Gaussian mixture (copula) models
#'
#' Easy and fast simulation of data from Gaussian mixture copula models (GMCM)
#' and Gaussian mixture models (GMM).
#'
#' @details
#' The functions provide simulation of \eqn{n} observations and
#' \eqn{d}-dimensional GMCMs and GMMs with provided parameters.
#' The \code{par} argument specifies the parameters of the Li et. al. (2011)
#' GMCM. The \code{theta} argument specifies an arbitrary GMCM of
#' Tewari et. al. (2011). Either one can be supplied. If both are missing,
#' random parameters are chosen for the general model.
#'
#' @aliases SimulateGMCMData SimulateGMMData
#' @param n A single integer giving the number of realizations (observations)
#'   drawn from the model. Default is 1000.
#' @param par A vector of parameters of length 4 where \code{par[1]} is the
#'   mixture proportion, \code{par[2]} is the mean, \code{par[3]} is the
#'   standard deviation, and \code{par[4]} is the correlation.
#' @param d The number of dimensions (or, equivalently, experiments) in the
#'   mixture distribution.
#' @param theta A list of parameters for the model as described in
#'   \code{\link{rtheta}}.
#' @param \dots In \code{SimulateGMCMData} the arguments are passed to
#'   \code{SimulateGMMData}.  In \code{SimulateGMMData} the arguments are passed
#'   to \code{\link{rtheta}}.
#' @return \code{SimulateGMCMData} returns a list of length 4 with elements:
#'   \item{u}{A matrix of the realized values of the GMCM.}
#'   \item{z}{A matrix of the latent GMM realizations.}
#'   \item{K}{An integer vector denoting the component from which the
#'     realization comes.}
#'   \item{theta}{A list containing the used parameters for the simulations
#'     with the format described in \code{\link{rtheta}}.}
#' @author Anders Ellern Bilgrau <anders.ellern.bilgrau@@gmail.com>
#' @seealso \code{\link{rtheta}}
#' @references
#'   Li, Q., Brown, J. B. J. B., Huang, H., & Bickel, P. J. (2011).
#'   Measuring reproducibility of high-throughput experiments. The Annals of
#'   Applied Statistics, 5(3), 1752-1779. doi:10.1214/11-AOAS466
#'
#'   Tewari, A., Giering, M. J., & Raghunathan, A. (2011). Parametric
#'   Characterization of Multimodal Distributions with Non-gaussian Modes.
#'   2011 IEEE 11th International Conference on Data Mining Workshops,
#'   286-292. doi:10.1109/ICDMW.2011.135
#' @examples
#' set.seed(2)
#'
#' # Simulation from the GMM
#' gmm.data1 <- SimulateGMMData(n = 200, m = 3, d = 2)
#' str(gmm.data1)
#'
#' # Plotting the simulated data
#' plot(gmm.data1$z, col = gmm.data1$K)
#'
#' # Simulation from the GMCM
#' gmcm.data1 <- SimulateGMCMData(n = 1000, m = 4, d = 2)
#' str(gmcm.data1)
#'
#' # Plotthe 2nd simulation
#' par(mfrow = c(1,2))
#' plot(gmcm.data1$z, col = gmcm.data1$K)
#' plot(gmcm.data1$u, col = gmcm.data1$K)
#'
#' # Simulation from the special case of GMCM
#' theta <- meta2full(c(0.7, 2, 1, 0.7), d = 3)
#' gmcm.data2 <- SimulateGMCMData(n = 5000, theta = theta)
#' str(gmcm.data2)
#'
#' # Plotting the 3rd simulation
#' par(mfrow=c(1,2))
#' plot(gmcm.data2$z, col = gmcm.data2$K)
#' plot(gmcm.data2$u, col = gmcm.data2$K)
#' @export
SimulateGMCMData <- function (n = 1000, par, d = 2, theta, ...) {
  if (missing(theta) & missing(par)) {
    theta <- rtheta(d = d, ...)
  } else if (missing(theta) & !missing(par)) {
    theta <- meta2full(par, d = d)
  }
  gmm.data <- SimulateGMMData(n = n, theta = theta, ...)
  u        <- pgmm.marginal(gmm.data$z, theta = theta)
  return(list(u = u, z = gmm.data$z, K = gmm.data$K, theta = theta))
}


