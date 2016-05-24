#' @rdname theta-utils
#' @description
#' \code{theta2tau} converts \eqn{\theta} to the transformation vector 
#' \eqn{\tau = (\mu_x, \sigma_x, \gamma, \delta, \alpha)}.
#' @seealso
#' \code{\link{beta2tau}}
#' @export
theta2tau <- function(theta = list(beta = c(0, 1)), distname, 
                      use.mean.variance = TRUE) {
  
  check_distname(distname)
  theta <- complete_theta(theta)
  tau <- beta2tau(beta = theta$beta, distname = distname,
                  use.mean.variance = use.mean.variance)

  if (length(theta$alpha) == 1) {
    tau["alpha"] <- theta$alpha
    tau <- tau[setdiff(names(tau), c("alpha_l", "alpha_r"))] 
  } else if (length(theta$alpha) == 2) {
    tau[c("alpha_l", "alpha_r")] <- theta$alpha    
    tau <- tau[setdiff(names(tau), c("alpha"))] 
  } else if (length(theta$alpha) > 2) {
    stop("alpha can't have more than 2 entries.")
  }
  
  if (length(theta$gamma) == 1) {
    tau["gamma"] <- theta$gamma
  }
  
  if (length(theta$delta) == 1) {
    tau["delta"] <- theta$delta
    tau <- tau[setdiff(names(tau), c("delta_l", "delta_r"))] 
  } else if (length(theta$delta) == 2) {
    tau[c("delta_l", "delta_r")] <- theta$delta
    tau <- tau[setdiff(names(tau), "delta")] 
  }
  return(tau)
}
