#' @rdname dpqr-zip
#' @name ZIP
#' @aliases dzip pzip qzip rzip
#' @title The Zero-Inflated Poisson Distribution
#' @description Density, distribution function, quantile function and random generation for the 
#' zero-inflated Poisson (ZIP) distribution with parameters \code{lambda} and \code{omega}.
#' @usage
#' dzip(x, lambda, omega, log = FALSE)
#' pzip(q, lambda, omega, lower.tail = TRUE, log.p = FALSE)
#' qzip(p, lambda, omega, lower.tail = TRUE, log.p = FALSE)
#' rzip(n, lambda, omega)
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param lambda vector of (non-negative) means.
#' @param omega zero-inflation parameter.
#' @param log,log.p logical; if TRUE, probabilities \code{p} are given as \code{log(p)}.
#' @param lower.tail logical; if TRUE (default), probabilities are \code{P[X <= x]}, otherwise, \code{P[X > x]}.
#' @return \code{dzip} gives the density, \code{pzip} gives the distribution function,
#' \code{qzip} gives the quantile function, and \code{rzip} generates random deviates.
#' @seealso \code{\link{dzinb}}, \code{\link{pzinb}}, \code{\link{qzinb}}, and \code{\link{rzinb}}
#' for the zero-inflated negative binomial (ZINB) distribution.
#' @examples 
#' dzip(x = 0:10, lambda = 1, omega = 0.5)
#' pzip(q = c(1, 5, 9), lambda = 1, omega = 0.5)
#' qzip(p = c(0.25, 0.50, 0.75), lambda = 1, omega = 0.5)
#' rzip(n = 100, lambda = 1, omega = 0.5)
#' @keywords distribution
#' @export dzip pzip qzip rzip
NULL

dzip <- function(x, lambda, omega, log = FALSE) {
  d <- omega * (x == 0) + (1 - omega) * dpois(x, lambda) 
  if(log == FALSE) {
    d
  } else if(log == TRUE) {
    log(d)
  }              
}

pzip <- function(q, lambda, omega, lower.tail = TRUE, log.p = FALSE) {
  p <- omega * (q >= 0) + (1 - omega) * ppois(q, lambda)
  if(lower.tail == FALSE) {
    p <- 1 - p
  }
  if(log.p == TRUE) {
    p <- log(p)
  }
  p
}

qzip <- function(p, lambda, omega, lower.tail = TRUE, log.p = FALSE) {
  if(lower.tail == FALSE) {
    p <- 1 - p
  }
  if(log.p == TRUE) {
    p <- exp(p)
  }
  qpois(pmax(0, (p - omega) / (1 - omega)), lambda)
}

rzip <- function(n, lambda, omega) {
  ifelse(rbinom(n, 1, omega), 0, rpois(n, lambda))
}


#' @rdname dpqr-zinb
#' @name ZINB
#' @aliases dzinb pzinb qzinb rzinb
#' @title The Zero-Inflated Negative Binomial Distribution
#' @description Density, distribution function, quantile function and random generation for the zero-inflated 
#' negative binomial (ZINB) distribution with parameters \code{k}, \code{lambda}, and \code{omega}.
#' @usage
#' dzinb(x, k, lambda, omega, log = FALSE)
#' pzinb(q, k, lambda, omega, lower.tail = TRUE, log.p = FALSE)
#' qzinb(p, k, lambda, omega, lower.tail = TRUE, log.p = FALSE)
#' rzinb(n, k, lambda, omega)
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param k dispersion parameter.
#' @param lambda vector of (non-negative) means.
#' @param omega zero-inflation parameter.
#' @param log,log.p logical; if TRUE, probabilities \code{p} are given as \code{log(p)}.
#' @param lower.tail logical; if TRUE (default), probabilities are \code{P[X <= x]}, otherwise, \code{P[X > x]}.
#' @return \code{dzinb} gives the density, \code{pzinb} gives the distribution function,
#' \code{qzinb} gives the quantile function, and \code{rzinb} generates random deviates.
#' @seealso \code{\link{dzip}}, \code{\link{pzip}}, \code{\link{qzip}}, and \code{\link{rzip}}
#' for the zero-inflated Poisson (ZIP) distribution.
#' @examples
#' dzinb(x = 0:10, k = 1, lambda = 1, omega = 0.5)
#' pzinb(q = c(1, 5, 9), k = 1, lambda = 1, omega = 0.5)
#' qzinb(p = c(0.25, 0.50, 0.75), k = 1, lambda = 1, omega = 0.5)
#' rzinb(n = 100, k = 1, lambda = 1, omega = 0.5)
#' @keywords distribution
#' @export dzinb pzinb qzinb rzinb
NULL

dzinb <- function(x, k, lambda, omega, log = FALSE) {
  d <- omega * (x == 0) + (1 - omega) * dnbinom(x, k, mu = lambda)
  if(log == FALSE) {
    d
  } else if(log == TRUE) {
    log(d)
  }
}

pzinb <- function(q, k, lambda, omega, lower.tail = TRUE, log.p = FALSE) {
  p <- omega * (q >= 0) + (1 - omega) * pnbinom(q, k, mu = lambda)
  if(lower.tail == FALSE) {
    p <- 1 - p
  }
  if(log.p == TRUE) {
    p <- log(p)
  }
  p
}

qzinb <- function(p, k, lambda, omega, lower.tail = TRUE, log.p = FALSE) {
  if(lower.tail == FALSE) {
    p <- 1 - p
  }
  if(log.p == TRUE) {
    p <- exp(p)
  }
  qnbinom(pmax(0, (p - omega) / (1 - omega)), k, mu = lambda)
}

rzinb <- function(n, k, lambda, omega) {
  ifelse(rbinom(n, 1, omega), 0, rnbinom(n, k, mu = lambda))
}
