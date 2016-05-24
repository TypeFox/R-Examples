#' @name rgw
#' @title Simulation of Generalized Waring values
#' @description Random generation of values from a Generalized Waring distribution with parameters \code{a}, \code{k} and \code{ro}.
#' @param n number of random values to return.
#' @param a vector of (non-negative) first parameters.
#' @param k vector of (non-negative) second parameters.
#' @param ro vector of (non-negative) third parameters.
#' @return \code{rgw} is an auxiliar function which generates random samples from a Generalized Waring distribution to be used in the simulated envelope called by \code{residuals}.
#' @importFrom stats rbeta rgamma rpois
#' @export

rgw <- function(n, a, k, ro){
  if (sum(a <= 0) | sum(k <= 0) | sum(ro <= 0)) stop('Parameters must be positive')
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol
  if ((any(!is.wholenumber(n)==TRUE) | any(n <= 0))) stop('n must be a positive integer')
  u <- rbeta(n, ro, k)
  v <- u / (1 - u)
  lam <- rgamma(n, a, v)
  salida <- rpois(n, lam)
  return(salida)
}
