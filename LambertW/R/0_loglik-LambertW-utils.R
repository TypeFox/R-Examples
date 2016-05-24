#' @title Log-Likelihood for Lambert W\eqn{\times} F RVs
#' @name loglik-LambertW-utils
#' 
#' @description
#' Evaluates the log-likelihood for \eqn{\theta} given observations \code{y}.
#' 
#' @details
#' For heavy-tail Lambert W\eqn{\times} F distributions (\code{type = "h"} or
#'     \code{type = "hh"}) the log-likelihood decomposes into an input
#'     log-likelihood plus a penalty term for transforming the data.
#' 
#' For skewed Lambert W \eqn{\times} F distributions this decomposition only
#'     exists for non-negative input RVs (e.g., \code{"exp"}onential,
#'     \code{"gamma"}, \code{"f"}, \ldots{}). If negative values are possible
#'     (\code{"normal"}, \code{"t"}, \code{"unif"}, \code{"cauchy"}, \ldots{})
#'     then \code{loglik_input} and \code{loglik_penalty} return \code{NA}, but
#'     the value of the output log-likelihood will still be returned correctly
#'     as \code{loglik.LambertW}.
#' 
#' See Goerg (2016) for details on the decomposition of the log-likelihood into
#'     a log-likelihood on the input parameters plus a penalty term for
#'     transforming the data.
#' 
#' @inheritParams common-arguments
#' @param y a numeric vector of real values (the observed data).
#' @return 
#' \code{loglik_input} and \code{loglik_penalty} return a scalar;
#' \code{loglik_LambertW} returns a list with 3 values:
#' 
#' \item{loglik.input}{ loglikelihood of \code{beta} given the transformed data,} 
#' \item{loglik.penalty}{ penalty for transforming the data,} 
#' \item{loglik.LambertW}{ total log-likelihood of \code{theta} given the observed data; 
#' if the former two values exist this is simply their sum.}
#' @keywords univar distribution
#' @examples
#' set.seed(1)
#' yy <- rLambertW(n = 1000, distname = "normal", 
#'                 theta = list(beta = c(0, 1), delta = 0.2))
#' loglik_penalty(tau = theta2tau(list(beta = c(1, 1), delta = c(0.2, 0.2)),
#'                                distname = "normal"), 
#'                y = yy, type = "hh")
#' # For a type = 's' Lambert W x F distribution with location family input
#' # such a decomposition doesn't exist; thus NA.
#' loglik_penalty(tau = theta2tau(list(beta = c(1, 1), gamma = 0.03), 
#'                                distname = "normal"),
#'                is.non.negative = FALSE,
#'                y = yy, type = "s") 
#' # For scale-family input it does exist
#' loglik_penalty(tau = theta2tau(list(beta = 1, gamma = 0.01), 
#'                                distname = "exp"),
#'                is.non.negative = TRUE,
#'                y = yy, type = "s") 
#'                
#' # evaluating the Gaussian log-likelihood
#' loglik_input(beta = c(0, 1), x = yy, distname = "normal") # built-in version
#' # or pass your own log pdf function
#' loglik_input(beta = c(0, 1), x = yy, distname = "user", 
#'              log.dX = function(xx, beta = beta) { 
#'                 dnorm(xx, mean = beta[1], sd = beta[2], log = TRUE)
#'              })
#' \dontrun{
#' # you must specify distname = 'user'; otherwise it does not work
#' loglik_input(beta = c(0, 1), x = yy, distname = "mydist", 
#'              log.dX = function(xx, beta = beta) { 
#'                 dnorm(xx, mean = beta[1], sd = beta[2], log = TRUE)
#'                 })
#' }
#' 
#' ### loglik_LambertW returns all three values
#' loglik_LambertW(theta = list(beta = c(1, 1), delta = c(0.09, 0.07)), 
#'                 y = yy, type = "hh", distname ="normal")
#' 
#' # can also take a flattend vector; must provide names though for delta
#' loglik_LambertW(theta = flatten_theta(list(beta = c(1, 1), 
#'                                           delta = c(delta_l = 0.09, 
#'                                                     delta_r = 0.07))), 
#'                 y = yy, type = "hh", distname ="normal")
#' 

NULL
