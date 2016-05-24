##' The Bridge Distribution
##'
##' Density, distribution function, quantile function and random generation for
##' the bridge distribution with parameter \code{scale}. See Wang and Louis (2003).
##'
##' If \code{scale} is omitted, the default
##' value \code{1/2} is assumed.
##'
##' The Bridge distribution parameterized by
##' \code{scale} has distribution function
##' \deqn{ }{F(q) = 1 - 1/(pi*scale) * (pi/2 - atan( (exp(scale*q) + cos(scale*pi)) / sin(scale*pi) ))}
##' and density
##' \deqn{ }{f(x) = 1/(2*pi) * sin(scale*pi) / (cosh(scale*x) + cos(scale*pi)).}
##'
##' The mean is \eqn{\mu}{0} and the variance is
##' \eqn{\pi^2 (\phi^{-2} - 1) / 3 }{pi^2 * (scale^{-2} - 1) / 3 }.
##'
##' @aliases Bridge bridge dbridge pbridge qbridge rbridge bridgedist
## @usage
## dbridge(x, scale = 1/2, log = FALSE) \cr
## pbridge(q, scale = 1/2, lower.tail = TRUE, log.p = FALSE) \cr
## qbridge(p, scale = 1/2, lower.tail = TRUE, log.p = FALSE) \cr
## rbridge(n, scale = 1/2)
##'
##' @param x,q vector of quantiles.
##' @param p vector of probabilities.
##' @param n number of observations. If \code{length(n) > 1}, the length is
##' taken to be the number required.
##' @param scale scale parameter. The scale must be between 0 and 1. A scale of 1/sqrt(1+3/pi^2) gives unit variance.
##' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
##' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X
##' \le x]}, otherwise, \eqn{P[X > x]}.
##' @return \code{dbridge} gives the density, \code{pbridge} gives the
##' distribution function, \code{qbridge} gives the quantile function, and
##' \code{rbridge} generates random deviates.
##'
##' The length of the result is determined by \code{n} for \code{rbridge}, and
##' is the maximum of the lengths of the numerical arguments for the other
##' functions.
##'
##' The numerical arguments other than \code{n} are recycled to the length of
##' the result.  Only the first elements of the logical arguments are used.
##' @note Consult the vignette for some figures comparing the normal, logistic, and bridge
##' distributions.
##'
##'
##' @seealso \link{Distributions} for other standard distributions.
##' @references Wang, Z. and Louis, T.A. (2003) Matching conditional and marginal shapes in binary random intercept models using
##' a bridge distribution function. \emph{Biometrika}, 90(4), 765-775. <DOI:10.1093/biomet/90.4.765>
##'
##' See also:
##'
##' Swihart, B.J., Caffo, B.S., and Crainiceanu, C.M. (2013). A Unifying Framework for Marginalized Random-Intercept Models of Correlated Binary Outcomes. \emph{International Statistical Review}, 82 (2), 275-295 1-22.  <DOI: 10.1111/insr.12035>
##'
##' Griswold, M.E., Swihart, B.J., Caffo, B.S and Zeger, S.L. (2013). Practical marginalized multilevel models. Stat, 2(1), 129-142. <DOI: 10.1002/sta4.22>
##'
##' Heagerty, P.J. (1999). Marginally specified logistic-normal models for longitudinal binary data. Biometrics, 55(3), 688-698.  <DOI: 10.1111/j.0006-341X.1999.00688.x>
##'
##' Heagerty, P.J. and Zeger, S.L. (2000). Marginalized multilevel models and likelihood inference (with comments and a rejoinder by the authors). Stat. Sci., 15(1), 1-26.  <DOI: 10.1214/ss/1009212671>
##'
##' @source \code{[dpq]bridge} are calculated directly from the definitions.
##'
##' \code{rbridge} uses inversion.
##' @keywords distribution
##' @examples
##'   ## Confirm unit variance for scale = 1/sqrt(1+3/pi^2)
##'   var(rbridge(1e5, scale = 1/sqrt(1+3/pi^2)))  # approximately 1
##'
##' @name Bridge
##'
##' @rdname Bridge
##' @export
dbridge <- function(x, scale = 1/2, log = FALSE){
  maxlength <- max(length(x), length(scale))
  x <- rep(x    , length=maxlength)
  scale <- rep(scale, length=maxlength)
  d=1/(2*pi) * sin(scale*pi) / (cosh(scale*x) + cos(scale*pi))
  if(log[1]) d = log(d)
  d
}
##' @rdname Bridge
##' @export
pbridge <- function(q, scale = 1/2, lower.tail = TRUE, log.p = FALSE){
  maxlength <- max(length(q), length(scale))
  q <- rep(q    , length=maxlength)
  scale <- rep(scale, length=maxlength)
  p=1 - 1/(pi*scale) * (pi/2 - atan( (exp(scale*q) + cos(scale*pi)) / sin(scale*pi) ))
  if(!lower.tail[1]) p = 1-p
  if(log.p[1]) p = log(p)
  p
}
##' @rdname Bridge
##' @export
qbridge <- function(p, scale = 1/2, lower.tail = TRUE, log.p = FALSE){
  maxlength <- max(length(p), length(scale))
  p <- rep(p    , length=maxlength)
  scale <- rep(scale, length=maxlength)
  if(log.p[1]) p = exp(p)
  if(!lower.tail[1]) p = 1-p
  q=1/scale * log( sin(scale*pi*p) / sin(scale*pi*(1-p)) )
  q
}
##' @rdname Bridge
##' @export
rbridge <- function(n, scale = 1/2){
  if(length(n) == 1){
    r <- qbridge(stats::runif(n), scale[1])
  }else{
    r <- qbridge(stats::runif(length(n)), rep(scale,length=length(n)))
  }
  r
}
