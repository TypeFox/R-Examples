#' Posterior class probabilities, local, and adjusted IDRs.
#'
#' Functions for computing posterior cluster probabilities (\code{get.prob})
#' in the general GMCM as well as local and
#' adjusted irreproducibility discovery rates (\code{get.IDR}) in the
#' special GMCM.
#'
#' @aliases get.IDR get.prob get.idr
#' @param x A \code{matrix} of observations where rows corresponds to features
#'   and columns to studies.
#' @param par A vector of length 4 where \code{par[1]} is mixture proportion of
#'   the irreproducible component, \code{par[2]} is the mean value,
#'   \code{par[3]} is the standard deviation, and \code{par[4]} is the
#'   correlation of the reproducible component.
#' @param threshold The threshold level of the IDR rate.
#' @param theta A list of parameters for the full model as described in
#'   \code{\link{rtheta}}.
#' @param \dots Arguments passed to \code{\link{qgmm.marginal}}.
#' @return
#' \code{get.IDR} returns a list of length 5 with elements:
#'   \item{idr}{A vector of the local idr values. I.e. the posterior
#'     probability that \code{x[i, ]} belongs to the irreproducible component.}
#'   \item{IDR}{A vector of the adjusted IDR values.}
#'   \item{l}{The number of reproducible features at the specified
#'     \code{threshold}.}
#'   \item{threshold}{The IDR threshold at which features are deemed
#'     reproducible.}
#'   \item{Khat}{A vector signifying whether the corresponding feature is
#'     reproducible or not.}
#' @author Anders Ellern Bilgrau <anders.ellern.bilgrau@@gmail.com>
#' @references
#'   Li, Q., Brown, J. B. J. B., Huang, H., & Bickel, P. J. (2011).
#'   Measuring reproducibility of high-throughput experiments. The Annals of
#'   Applied Statistics, 5(3), 1752-1779. doi:10.1214/11-AOAS466
#'
#'   Tewari, A., Giering, M., & Raghunathan, A. (2011). Parametric
#'   Characterization of Multimodal Distributions with Non-gaussian Modes. IEEE
#'   11th International Conference on Data Mining Workshops, 2011, 286-292.
#'   doi:10.1109/ICDMW.2011.135
#' @examples
#' set.seed(1123)
#'
#' # True parameters
#' true.par <- c(0.9, 2, 0.7, 0.6)
#'
#' # Simulation of data from the GMCM model
#' data <-  SimulateGMCMData(n = 1000, par = true.par, d = 2)
#'
#' # Initial parameters
#' init.par <- c(0.5, 1, 0.5, 0.9)
#'
#' # Nelder-Mead optimization
#' nm.par   <- fit.meta.GMCM(data$u, init.par = init.par, method = "NM")
#'
#' # Get IDR values
#' res <- get.IDR(data$u, nm.par, threshold = 0.05)
#'
#' # Plot results
#' plot(data$u, col = res$Khat, pch = c(3,16)[data$K])
#' @export
get.IDR <- function (x, par, threshold = 0.05, ...) {
  theta <- meta2full(par, d = ncol(x))
  idr  <- get.idr(x, theta, ...)
  ord  <- order(idr)
  IDR  <- cummean(idr[ord])
  if (any(IDR < threshold)) {
    l <- max(which(IDR < threshold))
  } else {
    l <- 0
  }
  IDR  <- IDR[order(ord)]
  Khat <- c(rep(2, l), rep(1, length(idr) - l))[order(ord)]
  return(list(idr = idr, IDR = IDR, l = l, threshold = threshold, Khat = Khat))
}
