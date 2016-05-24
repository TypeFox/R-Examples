#' @rdname theta-utils
#' 
#' @description
#' \code{get_theta_bounds} returns lower and upper bounds for \eqn{\theta} 
#' (necessary for optimization such as \code{\link{MLE_LambertW}}).
#' @param not.negative logical; if \code{TRUE} it sets the lower bounds for \code{alpha}
#' and \code{delta} to \code{0}.  Default: \code{FALSE}.
#' @return
#' \code{get_theta_bounds} returns a list containing two vectors:
#' \item{lower}{ flattened vector of lower bounds for valid \eqn{\theta}, }
#' \item{upper}{ flattened vector of upper bounds for valid \eqn{\theta}. }
#' 
#' @export
#' 
get_theta_bounds <- function(distname, beta, type = c("s", "h", "hh"),
                             not.negative = FALSE) {
  
  type <- match.arg(type)
  check_distname(distname)
  if (distname %in% c("normal", "cauchy")) {
    lb <- c(-Inf, 0)
    ub <- c(Inf, Inf)
  } else if (distname == "t") {
    lb <- c(-Inf, 0, 0)
    ub <- c(Inf, Inf, Inf)
  } else if (distname %in% c("exp", "chisq")) {
    lb <- c(0)
    ub <- c(Inf)
  } else if (distname %in% c("unif")) {
    lb <- c(-Inf, -Inf)
    ub <- c(Inf, Inf)
  } else {
    lb <- rep(-Inf, length(beta))
    ub <- rep(Inf, length(beta))
  }
  names(lb) <- names(ub) <- names(beta)
  
  if (type == "s") {
    lb <- c(lb, gamma = -Inf)
    ub <- c(ub, gamma = Inf)
  } else if (type == "h") {
    if (not.negative) {
      lb <- c(lb, alpha = 0, delta = 0)
    } else {
      lb <- c(lb, alpha = -Inf, delta = -Inf)
    }
    ub <- c(ub, alpha = Inf, delta = Inf)
  } else if (type == "hh") {
    if (not.negative) {
      lb <- c(lb, alpha_l = 0, alpha_r = 0, delta_l = 0, delta_r = 0)
    } else {
      lb <- c(lb, alpha_l = -Inf, alpha_r = -Inf, delta_l = -Inf, delta_r = -Inf)
    }

    ub <- c(ub, alpha_l = Inf, alpha_r = Inf, delta_l = Inf, delta_r = Inf)
  }
  
  out <- list(lower = lb, upper = ub)
  return(out)
} 
