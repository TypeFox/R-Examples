#' compute the dual constraint function: 1 if lambda <= l_prime, and a phi(c_1) otherwise
#'
#' @param  eta,gamma  prior mean and standard error (all inputs must be vectors of the same length)
#' @param  l_prime crossing points
#' @param  lambda the dual variable (a scalar)
#
#' @return Value of the dual constraint (vector of same length as the inputs)
#'
#' @export
#'
g_fun <- function(eta, gamma, lambda, l_prime) {

  J  <- length(eta)
  g_fun <- rep(0, J)

  for (i in 1:J) {
    if (lambda <= l_prime[i]) {
      g_fun[i] <- 1
    } else{
      g_fun[i] <- pnorm(Re(c_1(eta[i], gamma[i], lambda)))
    }
  }
  return(g_fun)
}
