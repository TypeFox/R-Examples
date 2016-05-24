#' @title Get bounds for gamma
#' @name get_gamma_bounds
#' 
#' @description
#' \code{get_gamma_bounds} returns lower and upper bounds for \eqn{\gamma}, so
#' that the observed data range falls within the theoretical bounds of the
#' support of the distribution. This is only important for location family
#' input.
#' 
#' @details Skewed Lambert W\eqn{\times} F distributions have
#'     parameter-dependent support for location family input.  Thus the
#'     parameter \eqn{\gamma} must be bounded such that the observed data is
#'     within the theoretical support of the distribution.  This theoretical
#'     bounds are determined by the Lambert W function (\code{\link{W}}), which
#'     has only real-valued solutions for \eqn{z \geq -1 / \exp(1)}.  Thus,
#'     \code{\link{W_gamma}} has real-valued solutions only for \eqn{z \geq -1 /
#'     \exp(1) \gamma} These lower and upper bounds are determined by minimum
#'     and maxiumum of the normalized data \eqn{\mathbf{z} = (\mathbf{y} -
#'     \mu_x) / \sigma_x}.
#' 
#' @inheritParams common-arguments
#' 
#' @return \code{get_gamma_bounds} returns a vector of length 2 with
#'     \code{"lower"} and \code{"upper"} bounds of \eqn{\gamma} given the range
#'     of \code{y}.
#' 
#' @export
#' 
get_gamma_bounds <- function(y, tau) {
  
  check_tau(tau)
  zz <- normalize_by_tau(y, tau)
  # if min(z) >= 0, then this is a not.negative random variable and we should
  # use c(0, Inf) as lower and upper bound
  if (min(zz) >= 0) { 
    bounds <- c("lower" = 0, "upper" = Inf)
  } else {
    bounds <- c("lower" = -1/exp(1)/max(zz), "upper" = -1/exp(1)/min(zz))
  }
  return(bounds)
}
