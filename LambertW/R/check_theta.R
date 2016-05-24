#' @rdname theta-utils
#' @description 
#' \code{check_theta} checks if \eqn{\theta = (\alpha, \boldsymbol \beta, \gamma, \delta)}
#' describes a well-defined Lambert W distribution.
#' @return
#' \code{check_theta} throws an error if list \code{theta} does not
#' define a proper Lambert W \eqn{\times} F distribution; 
#' does nothing otherwise.
#' @seealso
#' \code{\link{check_beta}}
#' @export
check_theta <- function(theta, distname)  {

  check_beta(theta$beta, distname = distname)
  
  if ("gamma" %in% names(theta)) {
    stopifnot(length(theta$gamma) == 1)
  } else { # add 'gamma' = 0 so checks below work
    theta$gamma <- 0
  }
  # check that alpha and delta are of same length; either 1 (for type 'h') or 
  # both (!) are vectors of length 2 (for type 'hh')
  if (!is.null(theta$delta) && !is.null(theta$alpha)) {
    if (length(theta$delta) != length(theta$alpha)) {
      print(paste("'delta':", theta$delta))
      print(paste("'alpha':", theta$alpha))
      stop("'delta' and 'alpha' have different length")
    }
  }
  
  
  if (theta$gamma != 0.0) {
    if (any(theta$delta != 0.0)) {
      stop(paste0("Both parameters are non-zero.",
                  "\n Either gamma = 0 (skewness) or delta = 0 (heavy-tails)."))
    } 
    if (!is.null(theta$alpha)) {
      if (any(theta$alpha != 1)) {
        stop("Alpha = ", unlist(theta$alpha), " and gamma = ", theta$gamma, ".",
             " For type ='s' skewed Lambert W distributions alpha == 1 or not present.")
      }
      if (any(theta$alpha <= 0)) {
        stop("'alpha' must be non-negative.")
      }

    }
  }
  }