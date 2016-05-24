#' find the crossing point of the two competing maxima
#'
#'@param eta mean
#'@param gamma marginal variance: gamma^2 = 1 + sigma^2, so gamma must be greater than 1
#'
#'@export
#'
find_crossing <- function(eta, gamma) {

  #source("c_1.R")
  f <- function(x) 1 - exp(x) - (pnorm((c_1(eta, gamma, x, 1) - eta) / gamma)
                                 - exp(x) * pnorm(c_1(eta, gamma, x, 1)))

  log_l <-  - log(gamma) - eta^2 / (2 * (gamma^2 - 1))

  if (f(log_l) > 0) {
    x0 <- c(log_l, 0) #initial interval
    x <- uniroot(f, x0)
    lambda <- exp(x$root)
  } else {
    #within numerical precision, lambda equals the starting point
    lambda <- exp(log_l)
  }
  return(lambda)
}
