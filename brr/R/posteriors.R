#' @name Posterior_lambda 
#' @rdname Posterior_lambda
#' @title Posterior distribution on the incidence rate in the treated group
#' @description Density and random generation for the posterior distribution on 
#' the rate in the treated group. The distribution function and the quantile function 
#' are not available.
#' @details The pdf of the posterior distribution of the incidence rate \eqn{\lambda} involves 
#' the Kummer confluent hypergeometric function of the second kind. 
#' 
#' 
#' @param lambda vector of quantiles 
#' @param a non-negative shape parameter of the Gamma prior distribution on \eqn{\mu}
#' @param c,d non-negative shape parameters of the prior distribution on \eqn{\phi} 
#' @param S sample size in treated group
#' @param x,y counts in the treated group and control group 
#' @param n number of observations to be simulated
#' @param ... other arguments passed to \code{\link{GIBDist}}
#' 
#' @return \code{dpost_lambda} gives the density, and \code{rpost_lambda} samples from 
#' the distribution, and \code{spost_lambda} gives a summary of the distribution.
#' 
#' @note \code{Posterior_lambda} is a generic name for the functions documented. 
#' 
#' @examples 
#' curve(dpost_lambda(x, 2, 2, 2, 20, 1, 10), from=0, to=0.4)
#' spost_lambda(2, 2, 2, 20, 1, 10)
NULL
#'
#' @rdname Posterior_lambda
#' @export 
dpost_lambda <- function(lambda, a, c, d, S, x, y, ...){
  a.post <- a+x+y
  c.post <- c+x
  d.post <- a+d+y
  return( dGIB(lambda, a.post, d.post, c.post, S, ...) )
}
#'
#' @rdname Posterior_lambda
#' @export 
rpost_lambda <- function(n, a, c, d, S, x, y){
  a.post <- a+x+y
  c.post <- c+x
  d.post <- a+d+y
  return( rGIB(n, a.post, d.post, c.post, S) )
}
#'
#' @rdname Posterior_lambda
#' @export 
spost_lambda <- function(a, c, d, S, x, y, ...){
  return( summary_GIB(a+x+y, a+d+y, c+x, S, ...) )
}

#' @name Posterior_mu 
#' @rdname Posterior_mu
#' @title Posterior distribution on the rate in the control group
#' @description Density and random generation for the posterior distribution on 
#' the rate in the control group. The distribution function and the quantile function 
#' are not available.
#' @details The pdf of the posterior distribution of the incidence rate \eqn{\mu} involves 
#' the Kummer confluent hypergeometric function of the second kind. 
#' 
#' @param mu vector of quantiles 
#' @param a,b non-negative shape and rate parameter of the Gamma prior distribution on \eqn{\mu}
#' @param c,d non-negative shape parameters of the prior distribution on \eqn{\phi} 
#' @param T sample size in control group 
#' @param x,y counts in the treated group and control group 
#' @param n number of observations to be simulated
#' @param ... other arguments passed to \code{\link{GIBDist}}
#' 
#' @return \code{dpost_mu} gives the density, \code{rpost_mu} samples from the 
#' distribution, and \code{spost_mu} gives a summary of the distribution.
#' 
#' @note \code{Posterior_mu} is a generic name for the functions documented. 
#' 
#' @examples 
#' curve(dpost_mu(x, 2, 2, 2, 2, 10, 3, 8), from=0, to=2)
#' spost_mu(2, 2, 2, 2, 10, 3, 8, output="pandoc")
#' 
NULL
#'
#' @rdname Posterior_mu
#' @export 
dpost_mu <- function(mu, a, b, c, d, T, x, y, ...){
  a.post <- a+x+y
  b.post <- b+T
  c.post <- c+x
  d.post <- a+d+y
  dGIB(mu, a.post, c.post, d.post, b.post, ...)
}
#'
#' @rdname Posterior_mu
#' @export 
rpost_mu <- function(n, a, b, c, d, T, x, y){
  a.post <- a+x+y
  b.post <- b+T
  c.post <- c+x
  d.post <- a+d+y
  return(rGIB(n, a.post, c.post, d.post, b.post))
}
#'
#' @rdname Posterior_mu
#' @export 
spost_mu <- function(a, b, c, d, T, x, y, ...){
  summary_GIB(a+x+y, c+x, a+d+y, T+b, ...)
}


#' @name Posterior_phi 
#' @rdname Posterior_phi
#' @title Posterior distribution on the relative risk and the vaccine efficacy
#' @description Density, distribution function, quantile function and random 
#' generation for the posterior distribution on relative risk or the vaccine efficacy.
#' @details The prior distribution on the relative risk \eqn{\phi} is the Beta2 distribution 
#' with shape parameters \eqn{c} and \eqn{d} and scale parameter \eqn{(T+b)/S}.
#' 
#' @param phi,VE,q vector of quantiles 
#' @param p vector of probabilities
#' @param a,b non-negative shape parameter and rate parameter of the prior Gamma distribution on the control incidence rate
#' @param c,d non-negative shape parameters of the prior distribution on the relative risk 
#' @param S,T sample sizes in control group and treated group
#' @param x,y counts in the treated group and control group 
#' @param n number of observations to be simulated
#' @param ... other arguments passed to \code{\link{Beta2Dist}}
#' 
#' @return \code{dpost_phi} gives the density, \code{ppost_phi} the distribution function, \code{qpost_phi} the quantile function, 
#' \code{rpost_phi} samples from the distribution, and \code{spost_phi} 
#' gives a summary of the distribution.
#' 
#' @note \code{Posterior_phi} is a generic name for the functions documented. 
#' 
#' @examples 
#' a <- 2; b <- 2; c <- 3; d <- 4; S <- 1; T <- 1; x <- 2; y <- 6
#' spost_phi(a, b, c, d, S, T, x, y, output="pandoc")
#' require(magrittr)
#' phi <- seq(0, 6, length.out=100)
#' phi %>% { plot(., dpost_phi(., a, b, c, d, S, T, x, y), type="l") }
#' phi %>% { lines(., dprior_phi(., b, c, d, S, T),  col="red") }
#' 
NULL
#'
#' @rdname Posterior_phi
#' @export 
dpost_phi <- function(phi, a, b, c, d, S, T, x, y, ...){
  c.post<-c+x ; d.post<-a+d+y
  dprior_phi(phi, b, c.post, d.post, S, T, ...)
}
#'
#' @rdname Posterior_phi
#' @export 
dpost_VE <- function(VE, a, b, c, d, S, T, x, y, ...){
  dpost_phi(1-VE, a, b, c, d, S, T, x, y, ...)
}
#'
#' @rdname Posterior_phi
#' @export 
ppost_phi <- function(q, a, b, c, d, S, T, x, y, ...){
  c.post<-c+x ; d.post<-a+d+y
  pprior_phi(q, b, c.post, d.post, S, T, ...)
}
#' @rdname Posterior_phi
#' @export 
ppost_VE <- function(q, a, b, c, d, S, T, x, y, ...){
  1 - ppost_phi(1-q, a, b, c, d, S, T, x, y, ...)
}
#'
#' @rdname Posterior_phi
#' @export 
qpost_phi <- function(p, a, b, c, d, S, T, x, y, ...){
  c.post <- c+x ; d.post <- a+d+y
  qprior_phi(p, b, c.post, d.post, S, T, ...)
}
#' @rdname Posterior_phi
#' @export 
qpost_VE <- function(p, a, b, c, d, S, T, x, y, ...){
  1 - qpost_phi(1-p, a, b, c, d, S, T, x, y, ...)
}
#'
#' @rdname Posterior_phi
#' @export 
rpost_phi <- function(n, a, b, c, d, S, T, x, y){
  c.post <- c+x ; d.post <- a+d+y
  rprior_phi(n, b, c.post, d.post, S, T)
}
#'
#' @rdname Posterior_phi
#' @export 
spost_phi <- function(a, b, c, d, S, T, x, y, ...){
  c.post <- c+x ; d.post <- a+d+y
  sprior_phi(b, c.post, d.post, S, T, ...)
}
