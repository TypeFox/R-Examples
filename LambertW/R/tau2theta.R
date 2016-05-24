#' @rdname tau-utils
#' @description
#' \code{tau2theta} converts \eqn{\tau} to the parameter list \eqn{\theta} 
#' (inverse of \code{\link{theta2tau}}).
#' @export
#' 

tau2theta <- function(tau, beta) {
  
  stopifnot(is.numeric(beta),
            is.numeric(tau))
  
  theta <- list(beta = beta, 
                alpha = 1, gamma = 0, delta = 0)
  
  if (!is.na(tau["alpha"])) {
    theta$alpha <- tau["alpha"]
  } 
  if (!is.na(tau["delta"])) {
    theta$delta <- tau["delta"]
  }
  if (any(!is.na(tau[c("delta_l", "delta_r")]))) {
    theta$delta <- c(tau["delta_l"], tau["delta_r"])
    names(theta$delta) <- c("delta_l", "delta_r")
  }
  if (!is.na(tau["gamma"])) {
    theta$gamma <- tau["gamma"]
  }
  theta <- complete_theta(theta)
  return(theta)
}