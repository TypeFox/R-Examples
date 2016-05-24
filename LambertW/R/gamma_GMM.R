#' @title Estimate gamma
#' @name gamma_GMM
#' 
#' @description
#' This function minimizes the Euclidean distance between the theoretical
#'     skewness of a skewed Lambert W x Gaussian random variable and the sample
#'     skewness of the back-transformed data \eqn{W_{\gamma}(\boldsymbol z)} as
#'     a function of \eqn{\gamma} (see References). Only an interative
#'     application of this function will give a good estimate of \eqn{\gamma}
#'     (see \code{\link{IGMM}}).
#' 
#' @param z a numeric vector of data values.
#' 
#' @param gamma.init starting value for \eqn{\gamma}; default:
#'     \code{\link{gamma_Taylor}}.
#' 
#' @param skewness.x theoretical skewness of the input \eqn{X}; default:
#'     \code{0}.
#'
#' @param robust logical; if \code{TRUE}, robust measure of asymmetry
#'     (\code{\link{medcouple_estimator}}) will be used; default: \code{FALSE}.
#'
#' @param tol a positive scalar; tolerance level for terminating the iterative
#'     algorithm; default: \code{.Machine$double.eps^0.25}.
#'
#' @param not.negative logical; if \code{TRUE}, the estimate for \eqn{\gamma} is
#'     restricted to non-negative reals, which is useful for scale-family
#'     Lambert W\eqn{\times} F random variables. Default: \code{FALSE}.
#'
#' @param optim.fct string; which R optimization function should be used.  By
#'     default it uses \code{\link[stats]{optimize}} which is about 8-10x faster
#'     than \code{\link[stats]{nlminb}}.
#'
#' @return 
#' A list with two elements: 
#' \item{gamma}{ scalar; optimal \eqn{\gamma}, } 
#' \item{iterations}{number of iterations (\code{NA} for \code{"optimize"}).}
#' @seealso 
#' \code{\link{delta_GMM}} for the heavy-tail version of this
#' function; \code{\link{medcouple_estimator}} for a robust measure of asymmetry;
#' \code{\link{IGMM}} for an iterative method to estimate all parameters
#' jointly.
#' @keywords optimize
#' @export
#' @examples
#' 
#' # highly skewed
#' y <- rLambertW(n = 1000, theta = list(beta = c(1, 2), gamma = 0.5), 
#'                distname = "normal") 
#' gamma_GMM(y, optim.fct = "nlminb")
#' gamma_GMM(y)
#' 
gamma_GMM <- function(z, skewness.x = 0, gamma.init = gamma_Taylor(z),
                      robust = FALSE, 
                      tol = .Machine$double.eps^0.25, not.negative = FALSE,
                      optim.fct = c("optimize", "nlminb")) {
  
  stopifnot(tol > 0,
            is.numeric(skewness.x),
            length(skewness.x) == 1)
  # optimize is about 8x faster than nlm(); but it does not show iterations
  optim.fct <- match.arg(optim.fct)
  
  .obj_fct <- function(gamma) {
    u.d <- W_gamma(z, gamma)
    if (anyNA(u.d)) {
      return(Inf)
    } else {
      if (!robust) {
        empirical.skewness <- skewness(u.d) 
      } else {
        empirical.skewness <- medcouple_estimator(u.d)
      }
      return(lp_norm(empirical.skewness - skewness.x, 2))
    }
  }

  bounds <- get_gamma_bounds(z, tau = c("mu_x" = 0, "sigma_x" = 1, gamma = 0))
  if (not.negative) {
    bounds[1] <- 0
  }
  lb <- bounds[1] + 1e-7
  ub <- bounds[2] - 1e-7
  if (ub == Inf) {
    ub <- 100
  }
    
  if (optim.fct == "nlminb") {
    fit <- nlminb(gamma.init, .obj_fct, lower = lb, upper = ub,
                  control = list(abs.tol = tol))
    out <- list(gamma = fit$par,
                iterations = fit$iterations)
  } else if (optim.fct == "optimize") {
    fit <- optimize(f = .obj_fct, interval = c(lb, ub), tol = tol)
    out <- list(gamma = fit$minimum,
                iterations = NA)
  }
  
  if (not.negative && lp_norm(out$gamma, p = 1) < 1e-6) {
    out$gamma <- round(out$gamma, 5) # to make almost zero to zero
  }
  return(out)
}
