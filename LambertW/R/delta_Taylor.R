#' @title Estimate of delta by Taylor approximation
#' 
#' @description
#' Computes an initial estimate of \eqn{\delta} based on the Taylor
#' approximation of the kurtosis of Lambert W \eqn{\times} Gaussian RVs. See
#' Details for the formula.
#' 
#' This is the initial estimate for \code{\link{IGMM}} and \code{\link{delta_GMM}}.
#' @details
#' The second order Taylor approximation of the theoretical kurtosis of a
#' heavy tail Lambert W x Gaussian RV around \eqn{\delta = 0} 
#' equals
#' 
#' \deqn{ \gamma_2(\delta) = 3 + 12 \delta + 66 \delta^2 + \mathcal{O}(\delta^3). }
#' 
#' Ignoring higher order terms, using the empirical estimate on the left hand side, and 
#' solving for \eqn{\delta} yields (positive root) 
#' \deqn{\widehat{\delta}_{Taylor} = \frac{1}{66} \cdot \left( \sqrt{66
#' \widehat{\gamma}_2(\mathbf{y}) - 162}-6 \right), } 
#' where \eqn{\widehat{\gamma}_2(\mathbf{y})} is the empirical kurtosis of \eqn{\mathbf{y}}.
#' 
#' Since the kurtosis is finite only for \eqn{\delta < 1/4},
#' \code{delta_Taylor} upper-bounds the returned estimate by \eqn{0.25}.
#' 
#' @param y a numeric vector of data values.
#' @param kurtosis.y kurtosis of \eqn{y}; default: empirical kurtosis of data \code{y}.
#' @param distname string; name of the distribution. Currently only supports \code{"normal"}.
#' @return 
#' scalar; estimated \eqn{\delta}.
#' 
#' @seealso \code{\link{IGMM}}  to estimate all parameters jointly.
#' @keywords optimize
#' @export
#' @examples
#' 
#' set.seed(2)
#' # a little heavy-tailed (kurtosis does exist)
#' y <- rLambertW(n = 1000, theta = list(beta = c(0, 1), delta = 0.2), 
#'                distname = "normal")
#' # good initial estimate since true delta=0.2 close to 0, and
#' # empirical kurtosis well-defined.
#' delta_Taylor(y) 
#' delta_GMM(y) # iterative estimate
#' 
#' y <- rLambertW(n = 1000, theta = list(beta = c(0, 1), delta = 1), 
#'                distname = "normal") # very heavy-tailed (like a Cauchy)
#' delta_Taylor(y) # bounded by 1/4 (as otherwise kurtosis does not exist)
#' delta_GMM(y) # iterative estimate
#' 
delta_Taylor <- function(y, kurtosis.y = kurtosis(y), distname = "normal") {
  stopifnot(is.numeric(kurtosis.y),
            length(kurtosis.y) == 1,
            kurtosis.y > 0)
  distname <- match.arg(distname)

  if (distname == "normal") {
    if (66 * kurtosis.y - 162 > 0) {
      delta.hat <- max(0, 1/66 * (sqrt(66 * kurtosis.y - 162) - 6))
      delta.hat <- min(delta.hat, 1/4)
    } else {
      delta.hat <- 0
    }
  } else {
    stop("Other distribution than 'normal' is not supported for the Taylor",
         "approximation.")
  }
  return(delta.hat)
} 
