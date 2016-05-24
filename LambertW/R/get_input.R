#' @title Back-transform Y to X
#' @name get_input
#' @aliases get.input
#' 
#' @description
#' \code{get_input} back-transforms the observed data \eqn{\boldsymbol y} to the
#'     (approximate) input data \eqn{\boldsymbol x_{\tau}} using the
#'     transformation vector \eqn{\tau = (\mu_x(\boldsymbol \beta),
#'     \sigma_x(\boldsymbol \beta), \gamma, \alpha, \delta)}.
#' 
#' Note that \code{get.input} should be deprecated; however, since it was
#'     explicitly referenced in Goerg (2011) I keep it here for future
#'     reference.  New code should use \code{get_input} exclusively.
#' 
#' @param y a numeric vector of data values or an object of class
#' \code{LambertW_fit}.
#' @param return.u should the normalized input be returned; default:
#' \code{FALSE}.
#' @param ... arguments passed to \code{get_input}.
#' @inheritParams common-arguments
#' @return 
#' The (approximated) input data vector \eqn{\widehat{\boldsymbol
#' x}_{\tau}}.
#' 
#' For \code{gamma != 0} it uses the principal branch solution
#' \code{\link{W_gamma}(z, branch = 0)} to get a unique input. 
#' 
#' For \code{gamma = 0} the back-transformation is bijective 
#' (for any \eqn{\delta \geq 0, \alpha \geq 0}).
#' 
#' If \code{return.u = TRUE}, then it returns a list with 2 vectors
#' \item{u}{centered and normalized input \eqn{\widehat{\boldsymbol u}_{\theta}},} 
#' \item{x}{input data \eqn{\widehat{\boldsymbol x}_{\theta}}.}
#' @keywords manip
#' @seealso
#' \code{\link{get_output}}
#' @export
#' @examples
#' 
#' set.seed(12)
#' # unskew very skewed data
#' y <- rLambertW(n = 1000, theta = list(beta = c(0, 1), gamma = 0.3), 
#'                distname = "normal")
#' test_normality(y)
#' fit.gmm <- IGMM(y, type="s")
#' 
#' x <- get_input(y, fit.gmm$tau)
#' # the same as
#' x <- get_input(fit.gmm)
#' test_normality(x) # symmetric Gaussian
#' 

get_input <- function(y, tau, return.u = FALSE) {
  if (inherits(y, "LambertW_fit")) {
    tau <- y$tau
    y <- y$data
  }
  tau <- complete_tau(tau)
  check_tau(tau)
  
  zz <- normalize_by_tau(y, tau)
  delta.values <- tau[grepl("delta", names(tau))]
  if (all(delta.values == 0) && tau["gamma"] == 0) {
      # no transformation
      uu <- zz
    } else if (tau["gamma"] != 0 && all(delta.values == 0)) {
      uu <- W_gamma(zz, gamma = tau["gamma"])
  } else if (tau["gamma"] == 0 && any(delta.values != 0)) {
    if ("delta_l" %in% names(delta.values)) {
      uu <- W_2delta_2alpha(zz, delta = tau[c("delta_l", "delta_r")], 
                            alpha = tau[c("alpha_l", "alpha_r")])
    } else if ("delta" %in% names(delta.values)) {
      uu <- W_delta_alpha(zz, delta = tau["delta"], alpha = tau['alpha'])
    }
  } else {
    stop("Only one of gamma or delta (or delta_l/delta_r) can be non-zero.")
  }
  xx <- normalize_by_tau(uu, tau, inverse = TRUE)  
  if (return.u) {
    out <- list(u = uu, x = xx)
    return(out)
  } else {
    return(xx)
  }
}

#' @rdname get_input
#' @export
get.input <- function(...) {
  warning("DEPRECATED: Please use get_input() instead of get.input().")
  return(get_input(...))
}
