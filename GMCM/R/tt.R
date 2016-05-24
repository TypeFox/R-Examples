#' Reparametrization of GMCM parameters
#'
#' These functions map the four GMCM parameters in the model of Li et. al.
#' (2011) and Tewari et. al. (2011) onto the real line and back. The mixture
#' proportion is logit transformed. The mean and standard deviation are log
#' transformed. The correlation is translated and scaled to the interval (0,1)
#' and logit transformed by \code{\link{rho.transform}}.
#'
#' The functions are used only in the wrapper to \code{optim} when the GMCM
#' log-likelihood is optimized.
#'
#' \code{par[1]} should be between 0 and 1. \code{par[2]} and \code{par[3]}
#' should be non-negative. If \code{positive.rho} is \code{FALSE},
#' \code{par[4]} should be between \eqn{-1/(d-1)}{-1/(d-1)} and 1. Otherwise,
#' \code{positive.rho} should be between 0 and 1.
#'
#' @aliases tt inv.tt
#' @param tpar A vector of length 4 of the transformed parameter values where
#'   \code{tpar[1]} corresponds to the mixture proportion, \code{tpar[2]} the
#'   mean, \code{tpar[3]} the standard deviation, and \code{tpar[4]} the
#'   correlation.
#' @param d The dimension of the space.
#' @param positive.rho is logical. If \code{TRUE}, the correlation is
#'   transformed by a simple \code{\link{logit}} transformation. If
#'   \code{FALSE} the
#' \code{\link{rho.transform}} is used.
#' @return A vector of the transformed or inversely transformed values of
#'   length 4.
#'
#'   \code{tt} returns \code{par} as described above.
#' @author Anders Ellern Bilgrau <anders.ellern.bilgrau@@gmail.com>
#' @references
#'   Li, Q., Brown, J. B. J. B., Huang, H., & Bickel, P. J. (2011).
#'   Measuring reproducibility of high-throughput experiments. The Annals of
#'   Applied Statistics, 5(3), 1752-1779. doi:10.1214/11-AOAS466
#'
#'   Tewari, A., Giering, M. J., & Raghunathan, A. (2011). Parametric
#'   Characterization of Multimodal Distributions with Non-gaussian Modes. 2011
#'   IEEE 11th International Conference on Data Mining Workshops, 286-292.
#'   doi:10.1109/ICDMW.2011.135
#' @examples
#' par <- c(pie1 = 0.3, mu = 2, sigma = 0.5, rho = 0.8)
#' tpar <- GMCM:::inv.tt(par, d = 3, positive.rho = FALSE)
#' GMCM:::tt(tpar, d = 3, positive.rho = FALSE)
#' @keywords internal
tt <- function(tpar, d, positive.rho) {
  par      <- NA
  par[1]   <- inv.logit(tpar[1])
  par[2:3] <- exp(tpar[2:3])
  if (positive.rho) {
    par[4] <- inv.logit(tpar[4])
  } else {
    par[4] <- inv.rho.transform(tpar[4], d)
  }
  return(par)
}
