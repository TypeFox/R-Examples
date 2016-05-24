#' @rdname theta-utils
#' @description
#' \code{theta2unbounded} transforms \eqn{\theta} from the bounded space to an 
#' unrestricted space (by \eqn{\log}-transformation on 
#' \eqn{\sigma_x}, \eqn{\delta}, and \eqn{\alpha}; note that this restricts
#' \eqn{\gamma \geq 0}, \eqn{\delta \geq 0}, and \eqn{\alpha \geq 0}.). 
#' 
#' @param inverse logical; if \code{TRUE}, it transforms the unbounded
#'     \code{theta} back to the original, bounded space. Default: \code{FALSE}.
#' 
#' @details Converting \code{theta} to an unbouded space is especially useful
#'     for optimization routines (like \code{\link[stats]{nlm}}), which can be
#'     performed over an unconstrained space. The obtained optimum can be
#'     converted back to the original space using the inverse transformation
#'     (set \code{inverse = TRUE} transforms it via \eqn{\exp}) -- this
#'     guarantees that the estimate satisfies non-negativity constraints (if
#'     required). The main advantage is that this avoids using optimization
#'     routines with boundary constraints -- since they are much slower compared
#'     to uncostrained optimization.
#' 
#' @export
theta2unbounded <- function(theta, distname, type = c("h", "hh", "s"), 
                            inverse = FALSE) {
  
  check_distname(distname)
  type <- match.arg(type)
  
  # check each element of theta, since
  
  if ('beta' %in% names(theta)) {
    if (distname %in% c("normal", "cauchy", "t")) {
      # only the scale parameter > 0
      if (inverse) {
        theta$beta[2] <- exp(theta$beta[2])
      } else {
        theta$beta[2] <- log(theta$beta[2])
      }
      
      if (distname == "t") {
        if (inverse) {
          theta$beta[3] <- exp(theta$beta[3])
        } else {
          theta$beta[3] <- log(theta$beta[3])
        }
      }
    } else if (distname %in% c("gamma", "f", "chisq", "exp", "beta")) {
      # all parameters > 0
      if (inverse) {
        theta$beta <- exp(theta$beta) 
      } else {
        theta$beta <- log(theta$beta)
      }    
    }
  }
  
  if ('delta' %in% names(theta)) {
    if (type %in% c("h", "hh")) {
      if (inverse) {
        theta$delta <- exp(theta$delta) 
      } else {
        theta$delta <- log(theta$delta)
      }
    }
  }
  
  if ("alpha" %in% names(theta)) {
    if (type %in% c("h", "hh")) {
      if (inverse) {
        theta$alpha <- exp(theta$alpha) 
      } else {
        theta$alpha <- log(theta$alpha)  
      }
    }
  }
  
  if ("gamma" %in% names(theta)) {
    if (!get_distname_family(distname)$location) {
      if (inverse) {
        theta$gamma <- exp(theta$gamma) 
      } else {
        theta$gamma <- log(theta$gamma)
      }
    }
  }
  return(theta)
} 
