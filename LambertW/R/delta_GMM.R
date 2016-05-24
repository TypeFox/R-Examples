#' @title Estimate delta
#' 
#' @description
#' 
#' This function minimizes the Euclidean distance between the sample kurtosis of
#'     the back-transformed data \eqn{W_{\delta}(\boldsymbol z)} and a
#'     user-specified target kurtosis as a function of \eqn{\delta} (see
#'     References).  Only an iterative application of this function will give a
#'     good estimate of \eqn{\delta} (see \code{\link{IGMM}}).
#' 
#' @param z a numeric vector of data values.
#' @inheritParams common-arguments
#' @param kurtosis.x theoretical kurtosis of the input X; default: \code{3}
#' (e.g., for \eqn{X \sim} Gaussian).
#' @param skewness.x theoretical skewness of the input X. Only used if \code{type = "hh"}; 
#' default: \code{0} (e.g., for \eqn{X \sim} symmetric).
#' @param delta.init starting value for optimization; default: \code{\link{delta_Taylor}}.
#' @param tol a positive scalar; tolerance level for terminating 
#' the iterative algorithm; default: \code{.Machine$double.eps^0.25}.
#' @param not.negative logical; if \code{TRUE} the estimate for \eqn{\delta} is
#' restricted to the non-negative reals. Default: \code{FALSE}.
#' @param optim.fct which R optimization function should be used. Either \code{'optimize'} 
#' (only for \code{type = 'h'} and if \code{not.negative = FALSE}) or \code{'nlm'}.  
#' Performance-wise there is no big difference.
#' @param lower,upper lower and upper bound for optimization if \code{optim.fct = 'optimize'}
#' and \code{not.negative = FALSE}. Default: \code{-1} and \code{3} 
#' (this covers most real-world heavy-tail scenarios).
#' @return 
#' A list with two elements: 
#' \item{delta}{ optimal \eqn{\delta} for data \eqn{z}, } 
#' \item{iterations}{number of iterations (\code{NA} for \code{'optimize'}).}
#' 
#' @seealso \code{\link{gamma_GMM}} for the skewed version of this function;
#' \code{\link{IGMM}} to estimate all parameters jointly.
#' @keywords optimize
#' @export
#' @examples
#' 
#' # very heavy-tailed (like a Cauchy)
#' y <- rLambertW(n = 1000, theta = list(beta = c(1, 2), delta = 1), 
#'                distname = "normal")
#' delta_GMM(y) # after the first iteration
#' 

delta_GMM <- function(z, type = c("h", "hh"),
                      kurtosis.x = 3, skewness.x = 0, 
                      delta.init = delta_Taylor(z), 
                      tol = .Machine$double.eps^0.25, 
                      not.negative = FALSE,
                      optim.fct = c("nlm", "optimize"),
                      lower = -1, upper = 3) {
  
  stopifnot(is.numeric(kurtosis.x),
            is.numeric(skewness.x),
            length(skewness.x) == 1,
            length(kurtosis.x) == 1,
            kurtosis.x > 0,
            length(delta.init) <= 2,
            tol > 0,
            lower < upper)
  
  optim.fct <- match.arg(optim.fct)
  type <- match.arg(type)

  if (type == "h") {
    .obj_fct <- function(delta) {
      if (not.negative) {
        # convert delta to > 0
        delta <- exp(delta)
      }
      u.g <- W_delta(z, delta = delta)

      if (anyNA(u.g) || any(is.infinite(u.g))) {
        return(lp_norm(kurtosis.x, 2))
      } else {
        empirical.kurtosis <- kurtosis(u.g)
        # for delta -> Inf, u.g can become (numerically) a constant vector
        # thus kurtosis(u.g) = NA.  In this case set empirical.kurtosis
        # to a very large value and continue.
        if (is.na(empirical.kurtosis)) {
          empirical.kurtosis <- 1e10
          warning("Kurtosis estimate was NA. ",
                  "Set to large value (", empirical.kurtosis, 
                  ") for optimization to continue.\n",
                  "Please double-check results (in particular the 'delta' ",
                  "estimate).")
        }
        return(lp_norm(empirical.kurtosis - kurtosis.x, 2))
      }
    }
  } else if (type == "hh") {
    .obj_fct <- function(delta) {
      if (not.negative) {
        # convert delta to > 0
        delta <- exp(delta)
      }
      u.g <- W_2delta(z, delta = delta)
      if (anyNA(u.g) || any(is.infinite(u.g))) {
        return(lp_norm(kurtosis.x, 2) + lp_norm(skewness.x * 2, 2))
      } else {
        empirical.kurtosis <- kurtosis(u.g)
        empirical.skewness <- skewness(u.g)
        
        # for delta -> Inf, u.g can become (numerically) a constant vector
        # thus kurtosis(u.g) = NA.  In this case set empirical.kurtosis
        # to a very large value and continue.
        if (is.na(empirical.kurtosis)) {
          empirical.kurtosis <- 1e10
          warning("Kurtosis estimate was NA. ",
                  "Set to large value (", empirical.kurtosis, 
                  ") for optimization to continue.\n",
                  "Please double-check results (in particular the 'delta' ",
                  "estimates).")
        }
        if (is.na(empirical.skewness)) {
          # make it skewed the same way as the input
          empirical.skewness <- 1e10 * (2 * as.numeric(skewness(z) > 0) - 1)
          warning("Skewness estimate was NA. ",
                  "Set to large value (", empirical.skewness, 
                  ") for optimization to continue.\n",
                  "Please double-check results (in particular the 'delta' ",
                  "estimates).")
        }
        return(lp_norm(empirical.kurtosis - kurtosis.x, 2) + 
                 lp_norm(empirical.skewness - skewness.x, 2))
      }
    }
    if (length(delta.init) == 1) {
      delta.init <- delta.init * c(1.1, 0.9)
    }
    if (skewness(z) > 0) {
      # revert lower and upper delta if skewness is positive
      delta.init <- rev(delta.init)
    }
  }
  if (not.negative) {
    # add 0.01 so delta.init = 0 (or c(0, 0)) can be a possible value for
    # initial value (delta_Taylor)
    delta.init <- log(delta.init + 0.001) 
  }
  
  out <- list()
  # use optimize only if non.negative is false
  if (length(delta.init) == 1 && optim.fct == "optimize" && !not.negative) {
      fit <- optimize(f = .obj_fct, interval = c(lower, upper), tol = tol)
      delta.hat <- fit$minimum
      out[["iterations"]] <- NA
  } else {
    fit <- nlm(f = .obj_fct, p = delta.init, print.level = 0, steptol = tol)
    delta.hat <- fit$estimate
    out[["iterations"]] <- fit$iterations
  }
  
  if (not.negative) {
    delta.hat <- exp(delta.hat)
    # round it to 6 digits, so that values like 1e-9 become 0
    if (lp_norm(delta.hat, 1) < 1e-7)
    delta.hat <- round(delta.hat, 6)
  }
  names(delta.hat) <- NULL
  out[["delta"]] <- delta.hat
  return(out)
} 
