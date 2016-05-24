#' Convert between parameter formats
#'
#' These functions converts the parameters between the general Gaussian
#' mixture (copula) model and the special GMCM.
#' Most functions of the GMCM packages use the \code{theta}
#' format described in \code{\link{rtheta}}.
#'
#' If a \code{theta} is supplied which is not on the form of Li et. al. (2011)
#' the output is coerced by simply picking the first elements of the first mean
#' vector and first covariance matrix as mean and standard deviation,
#' respectively.
#'
#' @aliases meta2full full2meta
#' @param theta A list of parameters for the full model. Formatted as described
#'   in \code{\link{rtheta}}.
#' @param par A vector of length 4 where \code{par[1]} is the probability of
#'   coming from the first component, \code{par[2]} is the mean value,
#'   \code{par[3]} is the standard deviation, and \code{par[4]} is the
#'   correlation of the reproducible component.
#' @return \code{full2meta} returns a numeric vector of length 4 formatted as
#'   \code{par}.
#' @author Anders Ellern Bilgrau <anders.ellern.bilgrau@@gmail.com>
#' @seealso \code{\link{rtheta}}
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
#' theta <- GMCM:::rtheta(m = 2, d = 2)
#' print(par <- full2meta(theta))
#' print(theta.special.case <- meta2full(par, d = 2))
#' @export
full2meta <- function(theta) {
  if (theta$m != 2) {
    stop("Too many components as m != 2. The special GMCM for meta-analysis",
         "only supports 2 components.")
  }
  par <- c("pie1"  = unname(theta$pie[1]),
           "mu"    = theta$mu[[2]][1],
           "sigma" = sqrt(theta$sigma[[2]][1,1]),
           "rho"   = cov2cor(theta$sigma[[2]])[1,2])
  if (par[4] <= -1/(theta$d  - 1))
    stop("correlation coefficient is not valid")
  return(par)
}
