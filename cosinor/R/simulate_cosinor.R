#' Simulate data from a cosinor model
#'
#' This function simulates data from a cosinor model with a single covariate,
#' where the time scale is month, and optionally
#' allows for single covariate effects on the mean,
#' amplitude, and acrophase.
#'
#' @param n Sample size
#' @param beta.mean Effect on the mean (intercept)
#' @param beta.amp Effect on the amplitude
#' @param beta.acro Effect on the acrophase
#'
#' @export
#'
simulate_cosinor <- function(n, beta.mean = 2, beta.amp = 0, beta.acro = 0){

  ttt <- runif(n, 0, 12)
  X <- rbinom(n, 1, .3)

  ## tranformations
  rrr <- cos(2 * pi * ttt / 12)
  sss <- sin(2 * pi * ttt / 12)

  ## generate data
  Y <- 30 + beta.mean * X + beta.amp * X * rrr + beta.acro * X * sss + rrr + 6 * sss + rnorm(n, sd = .1)
  data.frame(X = X, Y = Y, time = ttt)

}
