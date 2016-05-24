#' @rdname loglik-LambertW-utils
#' @description
#' \code{loglik_LambertW} computes the log-likelihood of \eqn{\theta}
#' for a Lambert W \eqn{\times} F distribution given observations \code{y}.
#' 
#' @param return.negative logical; if \code{TRUE} it returns the negative
#'     log-likelihood as a scalar (which is useful for numerical
#'     \emph{minimization} algorithms for \emph{maximum} likelihood estimation);
#'     otherwise it returns a list of input log-likelihood, penalty, and their
#'     sum = full likelihood. Default: \code{FALSE}.
#' 
#' @param flattened.theta.names vector of strings with names of flattened
#'     \code{theta}; this is necessary for optimization functions since they
#'     drop the names of a vector, but all functions in this package use names
#'     to select elements of (the flattened) \code{theta}.
#' 
#' @export
loglik_LambertW <- function(theta, y, distname, type, 
                            return.negative = FALSE,
                            flattened.theta.names = names(theta),
                            use.mean.variance = TRUE) {
  
  stopifnot(is.numeric(y))
  check_distname(distname)
  
  if (is.numeric(theta)) {
    names(theta) <- flattened.theta.names
    theta <- unflatten_theta(theta, distname = distname, type = type) 
  }
  tau <- theta2tau(theta, distname = distname, 
                   use.mean.variance = use.mean.variance)
  
  yy <- y
  xx <- get_input(yy, tau = tau)
  
  dist.family <- get_distname_family(distname)
  if (type == "s") {
    if (dist.family$is.non.negative) {
      out <- list(loglik.input = loglik_input(beta = theta$beta, x = xx, distname = distname),
                  loglik.penalty = loglik_penalty(tau = tau, y = yy, type = type,
                                                  is.non.negative = 
                                                      dist.family$is.non.negative))
      out$loglik.LambertW <- out$loglik.input + out$loglik.penalty
    } else {
      out <- list(loglik.input = NA,
                  loglik.penalty = NA,
                  loglik.LambertW = sum(dLambertW(yy, distname = distname, 
                                                  theta = theta, log = TRUE,
                                                  use.mean.variance = use.mean.variance)))
    }
  } else if (type %in% c("h", "hh")) {
    out <- list(loglik.input = loglik_input(beta = theta$beta, x = xx, distname = distname),
                loglik.penalty = loglik_penalty(tau = tau, y = yy, type = type))
    out$loglik.LambertW <- out$loglik.input + out$loglik.penalty
  }
  
  if (return.negative) {
    return(-out$loglik.LambertW)
  } else {
    return(out)
  }
}
