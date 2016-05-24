#' @name betagpd
#' 
#' @title Beta Bulk and GPD Tail Extreme Value Mixture Model
#'
#' @description Density, cumulative distribution function, quantile function and
#'   random number generation for the extreme value mixture model with beta for bulk
#'   distribution upto the threshold and conditional GPD above threshold. The parameters
#'   are the beta shape 1 \code{bshape1} and shape 2 \code{bshape2}, threshold \code{u}
#'   GPD scale \code{sigmau} and shape \code{xi} and tail fraction \code{phiu}.
#'
#' @param bshape1 beta shape 1 (positive)
#' @param bshape2 beta shape 2 (positive)
#' @param u       threshold over \eqn{(0, 1)}
#' @inheritParams normgpd
#' @inheritParams gpd
#' 
#' @details Extreme value mixture model combining beta distribution for the bulk
#' below the threshold and GPD for upper tail.
#' 
#' The user can pre-specify \code{phiu} 
#' permitting a parameterised value for the tail fraction \eqn{\phi_u}. Alternatively, when
#' \code{phiu=TRUE} the tail fraction is estimated as the tail fraction from the
#' beta bulk model.
#' 
#' The usual beta distribution is defined over \eqn{[0, 1]}, but this mixture is generally
#' not limited in the upper tail \eqn{[0,\infty]}, except for the usual upper tail 
#' limits for the GPD when \code{xi<0} discussed in \code{\link[evmix:gpd]{gpd}}. 
#' Therefore, the threshold is limited to \eqn{(0, 1)}.
#' 
#' The cumulative distribution function with tail fraction \eqn{\phi_u} defined by the
#' upper tail fraction of the beta bulk model (\code{phiu=TRUE}), upto the 
#' threshold \eqn{0 \le x \le u < 1}, given by:
#' \deqn{F(x) = H(x)}
#' and above the threshold \eqn{x > u}:
#' \deqn{F(x) = H(u) + [1 - H(u)] G(x)}
#' where \eqn{H(x)} and \eqn{G(X)} are the beta and conditional GPD
#' cumulative distribution functions (i.e. \code{pbeta(x, bshape1, bshape2)} and
#' \code{pgpd(x, u, sigmau, xi)}).
#' 
#' The cumulative distribution function for pre-specified \eqn{\phi_u}, upto the
#' threshold \eqn{0 \le x \le u < 1}, is given by:
#' \deqn{F(x) = (1 - \phi_u) H(x)/H(u)}
#' and above the threshold \eqn{x > u}:
#' \deqn{F(x) = \phi_u + [1 - \phi_u] G(x)}
#' Notice that these definitions are equivalent when \eqn{\phi_u = 1 - H(u)}.
#' 
#' See \code{\link[evmix:gpd]{gpd}} for details of GPD upper tail component and 
#' \code{\link[stats:Beta]{dbeta}} for details of beta bulk component.
#' 
#' @return \code{\link[evmix:betagpd]{dbetagpd}} gives the density, 
#' \code{\link[evmix:betagpd]{pbetagpd}} gives the cumulative distribution function,
#' \code{\link[evmix:betagpd]{qbetagpd}} gives the quantile function and 
#' \code{\link[evmix:betagpd]{rbetagpd}} gives a random sample.
#' 
#' @note All inputs are vectorised except \code{log} and \code{lower.tail}.
#' The main inputs (\code{x}, \code{p} or \code{q}) and parameters must be either
#' a scalar or a vector. If vectors are provided they must all be of the same length,
#' and the function will be evaluated for each element of vector. In the case of 
#' \code{\link[evmix:betagpd]{rbetagpd}} any input vector must be of length \code{n}.
#' 
#' Default values are provided for all inputs, except for the fundamentals 
#' \code{x}, \code{q} and \code{p}. The default sample size for 
#' \code{\link[evmix:betagpd]{rbetagpd}} is 1.
#' 
#' Missing (\code{NA}) and Not-a-Number (\code{NaN}) values in \code{x},
#' \code{p} and \code{q} are passed through as is and infinite values are set to
#' \code{NA}. None of these are not permitted for the parameters.
#' 
#' Error checking of the inputs (e.g. invalid probabilities) is carried out and
#' will either stop or give warning message as appropriate.
#' 
#' @references
#' \url{http://en.wikipedia.org/wiki/Beta_distribution}
#' 
#' \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
#' 
#' Scarrott, C.J. and MacDonald, A. (2012). A review of extreme value
#' threshold estimation and uncertainty quantification. REVSTAT - Statistical
#' Journal 10(1), 33-59. Available from \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}
#' 
#' MacDonald, A. (2012). Extreme value mixture modelling with medical and
#' industrial applications. PhD thesis, University of Canterbury, New Zealand.
#' \url{http://ir.canterbury.ac.nz/bitstream/10092/6679/1/thesis_fulltext.pdf}
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @seealso \code{\link[evmix:gpd]{gpd}} and \code{\link[stats:Beta]{dbeta}}
#' @aliases betagpd dbetagpd pbetagpd qbetagpd rbetagpd
#' @family  betagpd betagpdcon fbetagpd fbetagpdcon
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' par(mfrow = c(2, 2))
#' 
#' x = rbetagpd(1000, bshape1 = 1.5, bshape2 = 2, u = 0.7, phiu = 0.2)
#' xx = seq(-0.1, 2, 0.01)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-0.1, 2))
#' lines(xx, dbetagpd(xx, bshape1 = 1.5, bshape2 = 2, u = 0.7, phiu = 0.2))
#' 
#' # three tail behaviours
#' plot(xx, pbetagpd(xx, bshape1 = 1.5, bshape2 = 2, u = 0.7, phiu = 0.2), type = "l")
#' lines(xx, pbetagpd(xx, bshape1 = 1.5, bshape2 = 2, u = 0.7, phiu = 0.2, xi = 0.3), col = "red")
#' lines(xx, pbetagpd(xx, bshape1 = 1.5, bshape2 = 2, u = 0.7, phiu = 0.2, xi = -0.3), col = "blue")
#' legend("bottomright", paste("xi =",c(0, 0.3, -0.3)),
#'   col=c("black", "red", "blue"), lty = 1)
#' 
#' x = rbetagpd(1000, bshape1 = 2, bshape2 = 0.8, u = 0.7, phiu = 0.5)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-0.1, 2))
#' lines(xx, dbetagpd(xx, bshape1 = 2, bshape2 = 0.6, u = 0.7, phiu = 0.5))
#' 
#' plot(xx, dbetagpd(xx, bshape1 = 2, bshape2 = 0.8, u = 0.7, phiu = 0.5, xi=0), type = "l")
#' lines(xx, dbetagpd(xx, bshape1 = 2, bshape2 = 0.8, u = 0.7, phiu = 0.5, xi=-0.2), col = "red")
#' lines(xx, dbetagpd(xx, bshape1 = 2, bshape2 = 0.8, u = 0.7, phiu = 0.5, xi=0.2), col = "blue")
#' legend("topright", c("xi = 0", "xi = 0.2", "xi = -0.2"),
#'   col=c("black", "red", "blue"), lty = 1)
#' }
#' 
NULL

#' @export
#' @aliases betagpd dbetagpd pbetagpd qbetagpd rbetagpd
#' @rdname  betagpd

# probability density function for beta bulk with GPD for upper tail
dbetagpd <- function(x, bshape1 = 1, bshape2 = 1, u = qbeta(0.9, bshape1, bshape2),
  sigmau = sqrt(bshape1 * bshape2 / (bshape1 + bshape2)^2 / (bshape1 + bshape2 + 1)),
  xi = 0, phiu = TRUE, log = FALSE) {
    
  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.posparam(bshape1, allowvec = TRUE)
  check.posparam(bshape2, allowvec = TRUE)
  check.prob(u)
  check.posparam(sigmau, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.phiu(phiu, allowvec = TRUE)
  check.logic(log)

  n = check.inputn(c(length(x), length(bshape1), length(bshape2),
    length(u), length(sigmau), length(xi), length(phiu)), allowscalar = TRUE)

  if (any(is.infinite(x))) warning("infinite quantiles set to NA")

  x[is.infinite(x)] = NA # user will have to deal with infinite cases
  
  if ((min(u) <= 0) | (max(u) >= 1))
    stop("threshold must be between 0 and 1 (exclusive)")

  x = rep(x, length.out = n)
  bshape1 = rep(bshape1, length.out = n)
  bshape2 = rep(bshape2, length.out = n)
  u = rep(u, length.out = n)
  sigmau = rep(sigmau, length.out = n)
  xi = rep(xi, length.out = n)
  
  pu = pbeta(u, bshape1, bshape2)
  if (is.logical(phiu)) {
    phiu = 1 - pu
  } else {
    phiu = rep(phiu, length.out = n)
  }
  phib = (1 - phiu) / pu

  d = x # will pass through NA/NaN as entered
  
  whichb = which(x <= u)
  nb = length(whichb)
  whichu = which(x > u)
  nu = length(whichu)
  
  if (nb > 0) d[whichb] = log(phib[whichb]) + dbeta(x[whichb], bshape1[whichb], bshape2[whichb], log = TRUE)
  if (nu > 0) d[whichu] = log(phiu[whichu]) + dgpd(x[whichu], u[whichu], sigmau[whichu], xi[whichu], log = TRUE)

  if (!log) d = exp(d)

  d
}

#' @export
#' @aliases betagpd dbetagpd pbetagpd qbetagpd rbetagpd
#' @rdname  betagpd

# cumulative distribution function for beta bulk with GPD for upper tail
pbetagpd <- function(q, bshape1 = 1, bshape2 = 1, u = qbeta(0.9, bshape1, bshape2),
  sigmau = sqrt(bshape1 * bshape2 / (bshape1 + bshape2)^2 / (bshape1 + bshape2 + 1)),
  xi = 0, phiu = TRUE, lower.tail = TRUE) {

  # Check properties of inputs
  check.quant(q, allowna = TRUE, allowinf = TRUE)
  check.posparam(bshape1, allowvec = TRUE)
  check.posparam(bshape2, allowvec = TRUE)
  check.prob(u)
  check.posparam(sigmau, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.phiu(phiu, allowvec = TRUE)
  check.logic(lower.tail)

  n = check.inputn(c(length(q), length(bshape1), length(bshape2),
    length(u), length(sigmau), length(xi), length(phiu)), allowscalar = TRUE)

  if (any(is.infinite(q))) warning("infinite quantiles set to NA")

  q[is.infinite(q)] = NA # user will have to deal with infinite cases
  
  if ((min(u) <= 0) | (max(u) >= 1))
    stop("threshold must be between 0 and 1 (exclusive)")

  q = rep(q, length.out = n)
  bshape1 = rep(bshape1, length.out = n)
  bshape2 = rep(bshape2, length.out = n)
  u = rep(u, length.out = n)
  sigmau = rep(sigmau, length.out = n)
  xi = rep(xi, length.out = n)
  
  pu = pbeta(u, bshape1, bshape2)
  if (is.logical(phiu)) {
    phiu = 1 - pu
  } else {
    phiu = rep(phiu, length.out = n)
  }
  phib = (1 - phiu) / pu
    
  p = q # will pass through NA/NaN as entered
  
  whichb = which(q <= u)
  nb = length(whichb)
  whichu = which(q > u)
  nu = length(whichu)
  
  if (nb > 0) p[whichb] = phib[whichb] * pbeta(q[whichb], bshape1[whichb], bshape2[whichb])
  if (nu > 0) p[whichu] = 1 - phiu[whichu] + phiu[whichu] * pgpd(q[whichu], u[whichu], sigmau[whichu], xi[whichu])

  if (!lower.tail) p = 1 - p

  p
}

#' @export
#' @aliases betagpd dbetagpd pbetagpd qbetagpd rbetagpd
#' @rdname  betagpd

# inverse cumulative distribution function for beta bulk with GPD for upper tail
qbetagpd <- function(p, bshape1 = 1, bshape2 = 1, u = qbeta(0.9, bshape1, bshape2),
  sigmau = sqrt(bshape1 * bshape2 / (bshape1 + bshape2)^2 / (bshape1 + bshape2 + 1)),
  xi = 0, phiu = TRUE, lower.tail = TRUE) {

  # Check properties of inputs
  check.prob(p, allowna = TRUE)
  check.posparam(bshape1, allowvec = TRUE)
  check.posparam(bshape2, allowvec = TRUE)
  check.prob(u)
  check.posparam(sigmau, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.phiu(phiu, allowvec = TRUE)
  check.logic(lower.tail)

  n = check.inputn(c(length(p), length(bshape1), length(bshape2),
    length(u), length(sigmau), length(xi), length(phiu)), allowscalar = TRUE)
  
  if ((min(u) <= 0) | (max(u) >= 1))
    stop("threshold must be between 0 and 1 (exclusive)")

  if (!lower.tail) p = 1 - p

  p = rep(p, length.out = n)
  bshape1 = rep(bshape1, length.out = n)
  bshape2 = rep(bshape2, length.out = n)
  u = rep(u, length.out = n)
  sigmau = rep(sigmau, length.out = n)
  xi = rep(xi, length.out = n)
  
  pu = pbeta(u, bshape1, bshape2)
  if (is.logical(phiu)) {
    phiu = 1 - pu
  } else {
    phiu = rep(phiu, length.out = n)
  }
  phib = (1 - phiu) / pu
    
  q = p # will pass through NA/NaN as entered
  
  whichb = which(p <= (1 - phiu))
  nb = length(whichb)
  whichu = which(p > (1 - phiu))
  nu = length(whichu)

  if (nb > 0) q[whichb] = qbeta(p[whichb] / phib[whichb], bshape1[whichb], bshape2[whichb])
  if (nu > 0) q[whichu] = qgpd(p[whichu], u[whichu], sigmau[whichu], xi[whichu], phiu[whichu])

  q
}

#' @export
#' @aliases betagpd dbetagpd pbetagpd qbetagpd rbetagpd
#' @rdname  betagpd

# random number generation for beta bulk with GPD for upper tail
rbetagpd <- function(n = 1, bshape1 = 1, bshape2 = 1, u = qbeta(0.9, bshape1, bshape2),
  sigmau = sqrt(bshape1 * bshape2 / (bshape1 + bshape2)^2 / (bshape1 + bshape2 + 1)),
  xi = 0, phiu = TRUE) {

  # Check properties of inputs
  check.n(n)
  check.posparam(bshape1, allowvec = TRUE)
  check.posparam(bshape2, allowvec = TRUE)
  check.prob(u)
  check.posparam(sigmau, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.phiu(phiu, allowvec = TRUE)

  n = check.inputn(c(n, length(bshape1), length(bshape2), length(u), length(sigmau), length(xi), length(phiu)),
                   allowscalar = TRUE)

  if (any(xi == 1)) stop("shape cannot be 1")

  qbetagpd(runif(n), bshape1, bshape2, u, sigmau, xi, phiu)
}
