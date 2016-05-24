#' @title Utilities for the parameter vector of Lambert W\eqn{\times} F distributions
#' @name theta-utils
#' 
#' @description
#' These functions work with \eqn{\boldsymbol \theta = (\boldsymbol \beta, \gamma, \delta, \alpha)},
#' which fully parametrizes Lambert W\eqn{\times} F distributions.
#'
#' See Details for more background information on some functions.
#' @inheritParams common-arguments
#' @examples
#' 
#' \dontrun{
#' check_theta(theta = list(beta =  c(1, 1, -1)), distname = "t")
#' }
#' 
#' check_theta(theta = list(beta =  c(1, 1)), distname = "normal") # ok
#' 
#' params <- list(beta = c(2, 1), delta = 0.3) # alpha and gamma are missing
#' complete_theta(params) # added default values
#' 
#' params <- list(beta = c(2, 1), delta = 0.3, alpha = c(1, 2))
#' params <- complete_theta(params)
#' check_theta(params, distname = 'normal')
#' 
#' ###
#' x <- rnorm(1000)
#' get_initial_theta(x, distname = "normal", type = "h")
#' get_initial_theta(x, distname = "normal", type = "s")
#' 
#' # starting values for the skewed version of an exponential
#' y <- rLambertW(n = 1000, distname = "exp", beta = 2, gamma = 0.1)
#' get_initial_theta(y, distname = "exp", type = "s")
#' 
#' # starting values for the heavy-tailed version of a Normal = Tukey's h
#' y <- rLambertW(n = 1000, beta = c(2, 1), distname = "normal", delta = 0.2)
#' get_initial_theta(y, distname = "normal", type = "h")#' 
#' 
#' ###
#' get_theta_bounds(type = "hh", distname = "normal", beta = c(0, 1))
#' 
#' ### 
#' theta.restr <- theta2unbounded(list(beta = c(-1, 0.1), 
#'                                     delta = c(0.2, 0.2)), 
#'                                     distname = "normal")
#' theta.restr
#' # returns again the beta and delta from above
#' theta2unbounded(theta.restr, inverse = TRUE, distname = "normal") 
#' 
NULL
