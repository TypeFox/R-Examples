#' @name Prior_mu 
#' @rdname Prior_mu
#' @title Prior distribution on the rate in the control group
#' @description Density, distribution function, quantile function and random 
#' generation for the prior distribution on the rate in the control group.
#' @details The prior distribution on the rate \eqn{\mu} is the Gamma distribution 
#' with shape parameter \eqn{a} and rate parameter \eqn{b}
#' 
#' @param mu,q vector of quantiles 
#' @param p vector of probabilities
#' @param a,b non-negative shape parameter and rate parameter
#' @param n number of observations to be simulated
#' @param ... other arguments passed to \code{\link{GammaDist}} or to \code{\link{summary_gamma}}
#' 
#' @return \code{dprior_mu} gives the density, \code{pprior_mu} the distribution function, \code{qprior_mu} the quantile function, and \code{rprior_mu} samples from the distribution.
#' 
#' @note \code{Prior_mu} is a generic name for the functions documented. 
#' 
#' @examples 
#' curve(dprior_mu(x, 2, 2), from=0, to=3)
#' sprior_mu(2, 2, output="pandoc")
#' @importFrom stats dgamma pgamma qgamma rgamma
NULL
#'
#' @rdname Prior_mu
#' @export 
dprior_mu <- function(mu, a, b, ...){
  dgamma(mu, a, b, ...)
}
#'
#' @rdname Prior_mu
#' @export 
pprior_mu <- function(q, a, b, ...){
  pgamma(q, a, b, ...)
}
#'
#' @rdname Prior_mu
#' @export 
qprior_mu <- function(p, a, b, ...){
  qgamma(p, a, b, ...)
}
#'
#' @rdname Prior_mu
#' @export 
rprior_mu <- function(n, a, b, ...){
  rgamma(n, a, b, ...)
}
#'
#' @rdname Prior_mu
#' @export 
sprior_mu <- function(a, b, ...){
  summary_gamma(a, b, ...)
}

#' @name Prior_phi 
#' @rdname Prior_phi
#' @title Prior distribution on the relative risk and the vaccine efficacy
#' @description Density, distribution function, quantile function and random 
#' generation for the prior distribution on relative risk or the vaccine efficacy.
#' @details The prior distribution on the relative risk \eqn{\phi} is the Beta2 distribution 
#' with shape parameters \eqn{c} and \eqn{d} and scale parameter \eqn{(T+b)/S}.
#' 
#' @param phi,VE,q vector of quantiles 
#' @param p vector of probabilities
#' @param b non-negative rate parameter
#' @param c,d non-negative shape parameters 
#' @param S,T sample sizes in control group and treated group
#' @param n number of observations to be simulated
#' @param ... other arguments passed to \code{\link{Beta2Dist}}
#' 
#' @return \code{dprior_phi} gives the density, \code{pprior_phi} the distribution 
#' function, \code{qprior_phi} the quantile function, \code{rprior_phi} samples from 
#' the distribution, and \code{sprior_phi} gives a summary of the distribution.
#' 
#' @note \code{Prior_phi} is a generic name for the functions documented. 
#' 
#' @examples 
#' curve(dprior_phi(x, 2, 2, 2, 10, 10), from=0, to=7)
#' sprior_phi(2, 2, 2, 10, 10, output="pandoc")
#' 
NULL
#'
#' @rdname Prior_phi
#' @export 
dprior_phi<-function(phi, b, c, d, S, T, ...){
  scale <- (T+b)/S
  dbeta2(phi, c, d, scale, ...)
}
#' @rdname Prior_phi
#' @export 
dprior_VE<-function(VE, b, c, d, S, T, ...){
  dprior_phi(1-VE, b, c, d, S, T, ...)
}
#'
#' @rdname Prior_phi
#' @export 
pprior_phi<-function(q, b, c, d, S, T, ...){
  scale <- (T+b)/S
  pbeta2(q, c, d, scale, ...)
}
#' @rdname Prior_phi
#' @export 
pprior_VE<-function(q, b, c, d, S, T, ...){
  1-pprior_phi(1-q, b, c, d, S, T, ...)
}
#'
#' @rdname Prior_phi
#' @export 
qprior_phi <- function(p, b, c, d, S, T, ...){
  scale <- (T+b)/S
  qbeta2(p, c, d, scale, ...)
}
#' @rdname Prior_phi
#' @export 
qprior_VE <- function(p, b, c, d, S, T, ...){
  1-qprior_phi(1-p, b, c, d, S, T, ...)
}
#'
#' @rdname Prior_phi
#' @export 
rprior_phi <- function(n, b, c, d, S, T){
  rbeta2(n, c, d, scale=(T+b)/S)
}
#'
#' @rdname Prior_phi
#' @export 
sprior_phi <- function(b, c, d, S, T, ...){
  summary_beta2(c,d,scale=(T+b)/S, ...)
}


#' @name Prior_lambda 
#' @rdname Prior_lambda
#' @title Prior distribution on the incidence rate in the treated group
#' @description Density, distribution function (see Details) and random generation for the prior distribution on 
#' the rate in the treated group.  The prior distribution on the incidence rate \eqn{\lambda} is not to
#' be set by the user: it is induced by the user-specified prior on \eqn{\mu} 
#' and \eqn{\phi}.
#' @details The pdf of the prior distribution on the incidence rate \eqn{\lambda} involves 
#' the Kummer confluent hypergeometric function of the second kind. 
#' The cdf involves the generalized hypergeometric function. Its current implementation 
#' does not work when \code{a-c} is an integer.
#' 
#' @param lambda,q vector of quantiles 
#' @param a,b non-negative shape and rate parameter of the Gamma prior distribution on \eqn{\mu}
#' @param c,d non-negative shape parameters of the prior distribution on \eqn{\phi} 
#' @param S,T sample sizes in control group and treated group
#' @param n number of observations to be simulated
#' @param ... other arguments passed to \code{\link{genhypergeo}} through \code{\link{pGB2}}, 
#' such as \code{series=FALSE} to use the continued fraction expansion, 
#' or passed to \code{\link{summary_GB2}} (for \code{sprior_lambda})
#' 
#' @return \code{dprior_lambda} gives the density, \code{pprior_lambda} the distribution function 
#' (see Details), \code{rprior_lambda} samples from the distribution, and 
#' \code{rprior_lambda} gives a summary of the distribution.
#' 
#' @note \code{Prior_lambda} is a generic name for the functions documented. 
#' 
#' @examples 
#' curve(dprior_lambda(x, 2, 2, 2.5, 2, 10, 10), from=0, to=5)
#' sprior_lambda(2, 2, 2.5, 2, 10, 10)
NULL
#'
#' @rdname Prior_lambda
#' @importFrom gsl lnpoch lnbeta hyperg_U
#' @importFrom stats rgamma
#' @export 
dprior_lambda <- function(lambda, a, b, c, d, S, T){  
  return( dGB2(lambda, a=a, c=d, d=c, tau=b*S/(T+b)) )
}
#'
#' @rdname Prior_lambda
#' @export 
rprior_lambda <- function(n,a,b,c,d,S,T){
  return( rgamma(n,a,b) * rprior_phi(n, b, c, d, S, T) )
}
#'
#' @rdname Prior_lambda
#' @importFrom gsl lngamma lnpoch lnbeta hyperg_U
#' @importFrom hypergeo genhypergeo
#' @export 
pprior_lambda <- function(q, a, b, c, d, S, T, ...){ 
  if(a-c == floor(a-c)) stop("pprior_lambda() does not work when a-c is an integer.")
  return( pGB2(q, a=a, c=d, d=c, tau=b*S/(T+b), ...) )
}
#'
#' @rdname Prior_lambda
#' @export 
sprior_lambda <- function(a, b, c, d, S, T, ...){
  summary_GB2(a=a, c=d, d=c, tau=b*S/(T+b), ...)
}

