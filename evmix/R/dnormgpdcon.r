#' @name normgpdcon
#' 
#' @title Normal Bulk and GPD Tail Extreme Value Mixture Model with Single Continuity Constraint
#'
#' @description Density, cumulative distribution function, quantile function and
#'   random number generation for the extreme value mixture model with normal for bulk
#'   distribution upto the threshold and conditional GPD above threshold with continuity
#'   at threshold. The parameters
#'   are the normal mean \code{nmean} and standard deviation \code{nsd}, threshold \code{u}
#'   and GPD shape \code{xi} and tail fraction \code{phiu}.
#'
#' @inheritParams normgpd
#' 
#' @details Extreme value mixture model combining normal distribution for the bulk
#' below the threshold and GPD for upper tail with continuity at threshold.
#' 
#' The user can pre-specify \code{phiu} 
#' permitting a parameterised value for the tail fraction \eqn{\phi_u}. Alternatively, when
#' \code{phiu=TRUE} the tail fraction is estimated as the tail fraction from the
#' normal bulk model.
#' 
#' The cumulative distribution function with tail fraction \eqn{\phi_u} defined by the
#' upper tail fraction of the normal bulk model (\code{phiu=TRUE}), upto the 
#' threshold \eqn{x \le u}, given by:
#' \deqn{F(x) = H(x)}
#' and above the threshold \eqn{x > u}:
#' \deqn{F(x) = H(u) + [1 - H(u)] G(x)}
#' where \eqn{H(x)} and \eqn{G(X)} are the normal and conditional GPD
#' cumulative distribution functions (i.e. \code{pnorm(x, nmean, nsd)} and
#' \code{pgpd(x, u, sigmau, xi)}) respectively.
#' 
#' The cumulative distribution function for pre-specified \eqn{\phi_u}, upto the
#' threshold \eqn{x \le u}, is given by:
#' \deqn{F(x) = (1 - \phi_u) H(x)/H(u)}
#' and above the threshold \eqn{x > u}:
#' \deqn{F(x) = \phi_u + [1 - \phi_u] G(x)}
#' Notice that these definitions are equivalent when \eqn{\phi_u = 1 - H(u)}.
#' 
#' The continuity constraint means that \eqn{(1 - \phi_u) h(u)/H(u) = \phi_u g(u)}
#' where \eqn{h(x)} and \eqn{g(x)} are the normal and conditional GPD
#' density functions (i.e. \code{dnorm(x, nmean, nsd)} and
#' \code{dgpd(x, u, sigmau, xi)}) respectively. The resulting GPD scale parameter is then:
#' \deqn{\sigma_u = \phi_u H(u) / [1 - \phi_u] h(u)}.
#' In the special case of where the tail fraction is defined by the bulk model this reduces to
#' \deqn{\sigma_u = [1 - H(u)] / h(u)}. 
#' 
#' See \code{\link[evmix:gpd]{gpd}} for details of GPD upper tail component and 
#'\code{\link[stats:Normal]{dnorm}} for details of normal bulk component.
#' 
#' @return \code{\link[evmix:normgpdcon]{dnormgpdcon}} gives the density, 
#' \code{\link[evmix:normgpdcon]{pnormgpdcon}} gives the cumulative distribution function,
#' \code{\link[evmix:normgpdcon]{qnormgpdcon}} gives the quantile function and 
#' \code{\link[evmix:normgpdcon]{rnormgpdcon}} gives a random sample.
#' 
#' @note All inputs are vectorised except \code{log} and \code{lower.tail}.
#' The main inputs (\code{x}, \code{p} or \code{q}) and parameters must be either
#' a scalar or a vector. If vectors are provided they must all be of the same length,
#' and the function will be evaluated for each element of vector. In the case of 
#' \code{\link[evmix:normgpdcon]{rnormgpdcon}} any input vector must be of length \code{n}.
#' 
#' Default values are provided for all inputs, except for the fundamentals 
#' \code{x}, \code{q} and \code{p}. The default sample size for 
#' \code{\link[evmix:normgpdcon]{rnormgpdcon}} is 1.
#' 
#' Missing (\code{NA}) and Not-a-Number (\code{NaN}) values in \code{x},
#' \code{p} and \code{q} are passed through as is and infinite values are set to
#' \code{NA}. None of these are not permitted for the parameters.
#' 
#' Due to symmetry, the lower tail can be described by GPD by negating the quantiles. 
#' The normal mean \code{nmean} and GPD threshold \code{u} will also require negation.
#' 
#' Error checking of the inputs (e.g. invalid probabilities) is carried out and
#' will either stop or give warning message as appropriate.
#' 
#' @references
#' \url{http://en.wikipedia.org/wiki/Normal_distribution}
#' 
#' \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
#' 
#' Scarrott, C.J. and MacDonald, A. (2012). A review of extreme value
#' threshold estimation and uncertainty quantification. REVSTAT - Statistical
#' Journal 10(1), 33-59. Available from \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}
#' 
#' Behrens, C.N., Lopes, H.F. and Gamerman, D. (2004). Bayesian analysis of extreme
#' events with threshold estimation. Statistical Modelling. 4(3), 227-244.
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @seealso \code{\link[evmix:gpd]{gpd}} and \code{\link[stats:Normal]{dnorm}}
#' @aliases normgpdcon dnormgpdcon pnormgpdcon qnormgpdcon rnormgpdcon
#' @family  normgpd normgpdcon gng gngcon fnormgpd fnormgpdcon fgng fgngcon
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' par(mfrow = c(2, 2))
#' 
#' x = rnormgpdcon(1000)
#' xx = seq(-4, 6, 0.01)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-4, 6))
#' lines(xx, dnormgpdcon(xx))
#' 
#' # three tail behaviours
#' plot(xx, pnormgpdcon(xx), type = "l")
#' lines(xx, pnormgpdcon(xx, xi = 0.3), col = "red")
#' lines(xx, pnormgpdcon(xx, xi = -0.3), col = "blue")
#' legend("topleft", paste("xi =",c(0, 0.3, -0.3)),
#'   col=c("black", "red", "blue"), lty = 1)
#' 
#' x = rnormgpdcon(1000, phiu = 0.2)
#' xx = seq(-4, 6, 0.01)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-4, 6))
#' lines(xx, dnormgpdcon(xx, phiu = 0.2))
#' 
#' plot(xx, dnormgpdcon(xx, xi=0, phiu = 0.2), type = "l")
#' lines(xx, dnormgpdcon(xx, xi=-0.2, phiu = 0.2), col = "red")
#' lines(xx, dnormgpdcon(xx, xi=0.2, phiu = 0.2), col = "blue")
#' legend("topleft", c("xi = 0", "xi = 0.2", "xi = -0.2"),
#'   col=c("black", "red", "blue"), lty = 1)
#' }
#' 
NULL

#' @export
#' @aliases normgpdcon dnormgpdcon pnormgpdcon qnormgpdcon rnormgpdcon
#' @rdname  normgpdcon

# probability density function for normal bulk with GPD for upper tail
# with continuity at threshold
dnormgpdcon <- function(x, nmean = 0, nsd = 1, u = qnorm(0.9, nmean, nsd),
  xi = 0, phiu = TRUE, log = FALSE) {
  
  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(nmean, allowvec = TRUE)
  check.posparam(nsd, allowvec = TRUE)
  check.param(u, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.phiu(phiu, allowvec = TRUE)
  check.logic(log)

  n = check.inputn(c(length(x), length(nmean), length(nsd), length(u), length(xi), length(phiu)),
                   allowscalar = TRUE)

  if (any(is.infinite(x))) warning("infinite quantiles set to NA")

  x[is.infinite(x)] = NA # user will have to deal with infinite cases
  
  x = rep(x, length.out = n)
  nmean = rep(nmean, length.out = n)
  nsd = rep(nsd, length.out = n)
  u = rep(u, length.out = n)
  xi = rep(xi, length.out = n)
  
  pu = pnorm(u, nmean, nsd)
  if (is.logical(phiu)) {
    phiu = 1 - pu
  } else {
    phiu = rep(phiu, length.out = n)
  }
  phib = (1 - phiu) / pu
  
  sigmau = phiu / (phib * dnorm(u, nmean, nsd))
  
  check.posparam(sigmau, allowvec = TRUE)
  
  dnormgpd(x, nmean, nsd, u, sigmau, xi, phiu, log)
}

#' @export
#' @aliases normgpdcon dnormgpdcon pnormgpdcon qnormgpdcon rnormgpdcon
#' @rdname  normgpdcon

# cumulative distribution function for normal bulk with GPD for upper tail
# with continuity at threshold
pnormgpdcon <- function(q, nmean = 0, nsd = 1, u = qnorm(0.9, nmean, nsd), 
  xi = 0, phiu = TRUE, lower.tail = TRUE) {
  
  # Check properties of inputs
  check.quant(q, allowna = TRUE, allowinf = TRUE)
  check.param(nmean, allowvec = TRUE)
  check.posparam(nsd, allowvec = TRUE)
  check.param(u, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.phiu(phiu, allowvec = TRUE)
  check.logic(lower.tail)

  n = check.inputn(c(length(q), length(nmean), length(nsd), length(u), length(xi), length(phiu)),
                   allowscalar = TRUE)

  if (any(is.infinite(q))) warning("infinite quantiles set to NA")

  q[is.infinite(q)] = NA # user will have to deal with infinite cases
  
  q = rep(q, length.out = n)
  nmean = rep(nmean, length.out = n)
  nsd = rep(nsd, length.out = n)
  u = rep(u, length.out = n)
  xi = rep(xi, length.out = n)
  
  pu = pnorm(u, nmean, nsd)
  if (is.logical(phiu)) {
    phiu = 1 - pu
  } else {
    phiu = rep(phiu, length.out = n)
  }
  phib = (1 - phiu) / pu
  
  sigmau = phiu / (phib * dnorm(u, nmean, nsd))
  
  check.posparam(sigmau, allowvec = TRUE)

  pnormgpd(q, nmean, nsd, u, sigmau, xi, phiu, lower.tail)
}

#' @export
#' @aliases normgpdcon dnormgpdcon pnormgpdcon qnormgpdcon rnormgpdcon
#' @rdname  normgpdcon

# inverse cumulative distribution function for normal bulk with GPD for upper tail
# with continuity at threshold
qnormgpdcon <- function(p, nmean = 0, nsd = 1, u = qnorm(0.9, nmean, nsd), 
  xi = 0, phiu = TRUE, lower.tail = TRUE) {
  
  # Check properties of inputs
  check.prob(p, allowna = TRUE)
  check.param(nmean, allowvec = TRUE)
  check.posparam(nsd, allowvec = TRUE)
  check.param(u, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.phiu(phiu, allowvec = TRUE)
  check.logic(lower.tail)

  n = check.inputn(c(length(p), length(nmean), length(nsd), length(u), length(xi), length(phiu)),
                   allowscalar = TRUE)
    
  p = rep(p, length.out = n)
  nmean = rep(nmean, length.out = n)
  nsd = rep(nsd, length.out = n)
  u = rep(u, length.out = n)
  xi = rep(xi, length.out = n)
  
  pu = pnorm(u, nmean, nsd)
  if (is.logical(phiu)) {
    phiu = 1 - pu
  } else {
    phiu = rep(phiu, length.out = n)
  }
  phib = (1 - phiu) / pu
    
  sigmau = phiu / (phib * dnorm(u, nmean, nsd))
  
  check.posparam(sigmau, allowvec = TRUE)
    
  qnormgpd(p, nmean, nsd, u, sigmau, xi, phiu, lower.tail)
}

#' @export
#' @aliases normgpdcon dnormgpdcon pnormgpdcon qnormgpdcon rnormgpdcon
#' @rdname  normgpdcon

# random number generation for normal bulk with GPD for upper tail
# with continuity at threshold
rnormgpdcon <- function(n = 1, nmean = 0, nsd = 1, u = qnorm(0.9, nmean, nsd),
  xi = 0, phiu = TRUE) {
  
  # Check properties of inputs
  check.n(n)
  check.param(nmean, allowvec = TRUE)
  check.posparam(nsd, allowvec = TRUE)
  check.param(u, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.phiu(phiu, allowvec = TRUE)

  n = check.inputn(c(n, length(nmean), length(nsd), length(u), length(xi), length(phiu)), allowscalar = TRUE)
  
  if (any(xi == 1)) stop("shape cannot be 1")
  
  qnormgpdcon(runif(n), nmean, nsd, u, xi, phiu)
}
                                                                              