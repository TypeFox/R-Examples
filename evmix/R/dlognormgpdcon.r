#' @name lognormgpdcon
#' 
#' @title Log-Normal Bulk and GPD Tail Extreme Value Mixture Model with Single Continuity Constraint
#'
#' @description Density, cumulative distribution function, quantile function and
#'   random number generation for the extreme value mixture model with log-normal for bulk
#'   distribution upto the threshold and conditional GPD above threshold with continuity
#'   at threshold. The parameters
#'   are the log-normal mean \code{lnmean} and standard deviation \code{lnsd}, threshold \code{u}
#'   GPD shape \code{xi} and tail fraction \code{phiu}.
#'
#' @inheritParams lognormgpd
#' 
#' @details Extreme value mixture model combining log-normal distribution for the bulk
#' below the threshold and GPD for upper tailwith continuity
#'   at threshold.
#' 
#' The user can pre-specify \code{phiu} 
#' permitting a parameterised value for the tail fraction \eqn{\phi_u}. Alternatively, when
#' \code{phiu=TRUE} the tail fraction is estimated as the tail fraction from the
#' log-normal bulk model.
#' 
#' The cumulative distribution function with tail fraction \eqn{\phi_u} defined by the
#' upper tail fraction of the log-normal bulk model (\code{phiu=TRUE}), upto the 
#' threshold \eqn{0 < x \le u}, given by:
#' \deqn{F(x) = H(x)}
#' and above the threshold \eqn{x > u}:
#' \deqn{F(x) = H(u) + [1 - H(u)] G(x)}
#' where \eqn{H(x)} and \eqn{G(X)} are the log-normal and conditional GPD
#' cumulative distribution functions (i.e. \code{plnorm(x, lnmean, lnsd)} and
#' \code{pgpd(x, u, sigmau, xi)}) respectively.
#' 
#' The cumulative distribution function for pre-specified \eqn{\phi_u}, upto the
#' threshold \eqn{0 < x \le u}, is given by:
#' \deqn{F(x) = (1 - \phi_u) H(x)/H(u)}
#' and above the threshold \eqn{x > u}:
#' \deqn{F(x) = \phi_u + [1 - \phi_u] G(x)}
#' Notice that these definitions are equivalent when \eqn{\phi_u = 1 - H(u)}.
#' 
#' The log-normal is defined on the positive reals, so the threshold must be positive.
#' 
#' The continuity constraint means that \eqn{(1 - \phi_u) h(u)/H(u) = \phi_u g(u)}
#' where \eqn{h(x)} and \eqn{g(x)} are the log-normal and conditional GPD
#' density functions (i.e. \code{dlnorm(x, lnmean, lnsd)} and
#' \code{dgpd(x, u, sigmau, xi)}) respectively. The resulting GPD scale parameter is then:
#' \deqn{\sigma_u = \phi_u H(u) / [1 - \phi_u] h(u)}.
#' In the special case of where the tail fraction is defined by the bulk model this reduces to
#' \deqn{\sigma_u = [1 - H(u)] / h(u)}. 
#' 
#' See \code{\link[evmix:gpd]{gpd}} for details of GPD upper tail component and 
#'\code{\link[stats:Lognormal]{dlnorm}} for details of log-normal bulk component.
#' 
#' @return \code{\link[evmix:lognormgpdcon]{dlognormgpdcon}} gives the density, 
#' \code{\link[evmix:lognormgpdcon]{plognormgpdcon}} gives the cumulative distribution function,
#' \code{\link[evmix:lognormgpdcon]{qlognormgpdcon}} gives the quantile function and 
#' \code{\link[evmix:lognormgpdcon]{rlognormgpdcon}} gives a random sample.
#' 
#' @note All inputs are vectorised except \code{log} and \code{lower.tail}.
#' The main inputs (\code{x}, \code{p} or \code{q}) and parameters must be either
#' a scalar or a vector. If vectors are provided they must all be of the same length,
#' and the function will be evaluated for each element of vector. In the case of 
#' \code{\link[evmix:lognormgpdcon]{rlognormgpdcon}} any input vector must be of length \code{n}.
#' 
#' Default values are provided for all inputs, except for the fundamentals 
#' \code{x}, \code{q} and \code{p}. The default sample size for 
#' \code{\link[evmix:lognormgpdcon]{rlognormgpdcon}} is 1.
#' 
#' Missing (\code{NA}) and Not-a-Number (\code{NaN}) values in \code{x},
#' \code{p} and \code{q} are passed through as is and infinite values are set to
#' \code{NA}. None of these are not permitted for the parameters.
#' 
#' Error checking of the inputs (e.g. invalid probabilities) is carried out and
#' will either stop or give warning message as appropriate.
#' 
#' @references
#' \url{http://en.wikipedia.org/wiki/Log-normal_distribution}
#' 
#' \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
#' 
#' Scarrott, C.J. and MacDonald, A. (2012). A review of extreme value
#' threshold estimation and uncertainty quantification. REVSTAT - Statistical
#' Journal 10(1), 33-59. Available from \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}
#' 
#' Solari, S. and Losada, M.A. (2004). A unified statistical model for
#' hydrological variables including the selection of threshold for the peak over
#' threshold method. Water Resources Research. 48, W10541.
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @seealso \code{\link[evmix:gpd]{gpd}} and \code{\link[stats:Lognormal]{dlnorm}}
#' @aliases lognormgpdcon dlognormgpdcon plognormgpdcon qlognormgpdcon rlognormgpdcon
#' @family  lognormgpd lognormgpdcon normgpd normgpdcon flognormgpd flognormgpdcon fnormgpd fnormgpdcon
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' par(mfrow = c(2, 2))
#' 
#' x = rlognormgpdcon(1000)
#' xx = seq(-1, 10, 0.01)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-1, 10))
#' lines(xx, dlognormgpdcon(xx))
#' 
#' # three tail behaviours
#' plot(xx, plognormgpdcon(xx), type = "l")
#' lines(xx, plognormgpdcon(xx, xi = 0.3), col = "red")
#' lines(xx, plognormgpdcon(xx, xi = -0.3), col = "blue")
#' legend("bottomright", paste("xi =",c(0, 0.3, -0.3)),
#'   col=c("black", "red", "blue"), lty = 1)
#' 
#' x = rlognormgpdcon(1000, u = 2, phiu = 0.2)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-1, 10))
#' lines(xx, dlognormgpdcon(xx, u = 2, phiu = 0.2))
#' 
#' plot(xx, dlognormgpdcon(xx, u = 2, xi=0, phiu = 0.2), type = "l")
#' lines(xx, dlognormgpdcon(xx, u = 2, xi=-0.2, phiu = 0.2), col = "red")
#' lines(xx, dlognormgpdcon(xx, u = 2, xi=0.2, phiu = 0.2), col = "blue")
#' legend("topright", c("xi = 0", "xi = 0.2", "xi = -0.2"),
#'   col=c("black", "red", "blue"), lty = 1)
#' }
#' 
NULL

#' @export
#' @aliases lognormgpdcon dlognormgpdcon plognormgpdcon qlognormgpdcon rlognormgpdcon
#' @rdname  lognormgpdcon

# probability density function for log-normal bulk with GPD for upper tail
# with continuity at threshold
dlognormgpdcon <- function(x, lnmean = 0, lnsd = 1, u = qlnorm(0.9, lnmean, lnsd),
  xi = 0, phiu = TRUE, log = FALSE) {

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(lnmean, allowvec = TRUE)
  check.posparam(lnsd, allowvec = TRUE)
  check.posparam(u, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.phiu(phiu, allowvec = TRUE)
  check.logic(log)

  n = check.inputn(c(length(x), length(lnmean), length(lnsd), length(u), length(xi), length(phiu)),
                   allowscalar = TRUE)

  if (any(is.infinite(x))) warning("infinite quantiles set to NA")

  x[is.infinite(x)] = NA # user will have to deal with infinite cases
 
  x = rep(x, length.out = n)
  lnmean = rep(lnmean, length.out = n)
  lnsd = rep(lnsd, length.out = n)
  u = rep(u, length.out = n)
  xi = rep(xi, length.out = n)
  
  pu = plnorm(u, lnmean, lnsd)
  if (is.logical(phiu)) {
    phiu = 1 - pu
  } else {
    phiu = rep(phiu, length.out = n)
  }
  phib = (1 - phiu) / pu
  
  sigmau = phiu / (phib * dlnorm(u, lnmean, lnsd))
  
  check.posparam(sigmau, allowvec = TRUE)
  
  dlognormgpd(x, lnmean, lnsd, u, sigmau, xi, phiu, log)
}

#' @export
#' @aliases lognormgpdcon dlognormgpdcon plognormgpdcon qlognormgpdcon rlognormgpdcon
#' @rdname  lognormgpdcon

# cumulative distribution function for log-normal bulk with GPD for upper tail
# with continuity at threshold
plognormgpdcon <- function(q, lnmean = 0, lnsd = 1, u = qlnorm(0.9, lnmean, lnsd), 
  xi = 0, phiu = TRUE, lower.tail = TRUE) {

  # Check properties of inputs
  check.quant(q, allowna = TRUE, allowinf = TRUE)
  check.param(lnmean, allowvec = TRUE)
  check.posparam(lnsd, allowvec = TRUE)
  check.posparam(u, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.phiu(phiu, allowvec = TRUE)
  check.logic(lower.tail)

  n = check.inputn(c(length(q), length(lnmean), length(lnsd), length(u), length(xi), length(phiu)),
                   allowscalar = TRUE)

  if (any(is.infinite(q))) warning("infinite quantiles set to NA")

  q[is.infinite(q)] = NA # user will have to deal with infinite cases

  q = rep(q, length.out = n)
  lnmean = rep(lnmean, length.out = n)
  lnsd = rep(lnsd, length.out = n)
  u = rep(u, length.out = n)
  xi = rep(xi, length.out = n)
  
  pu = plnorm(u, lnmean, lnsd)
  if (is.logical(phiu)) {
    phiu = 1 - pu
  } else {
    phiu = rep(phiu, length.out = n)
  }
  phib = (1 - phiu) / pu
  
  sigmau = phiu / (phib * dlnorm(u, lnmean, lnsd))
  
  check.posparam(sigmau, allowvec = TRUE)

  plognormgpd(q, lnmean, lnsd, u, sigmau, xi, phiu, lower.tail)
}

#' @export
#' @aliases lognormgpdcon dlognormgpdcon plognormgpdcon qlognormgpdcon rlognormgpdcon
#' @rdname  lognormgpdcon

# inverse cumulative distribution function for log-normal bulk with GPD for upper tail
# with continuity at threshold
qlognormgpdcon <- function(p, lnmean = 0, lnsd = 1, u = qlnorm(0.9, lnmean, lnsd),
  xi = 0, phiu = TRUE, lower.tail = TRUE) {

  # Check properties of inputs
  check.prob(p, allowna = TRUE)
  check.param(lnmean, allowvec = TRUE)
  check.posparam(lnsd, allowvec = TRUE)
  check.posparam(u, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.phiu(phiu, allowvec = TRUE)
  check.logic(lower.tail)

  n = check.inputn(c(length(p), length(lnmean), length(lnsd), length(u), length(xi), length(phiu)),
                   allowscalar = TRUE)

  p = rep(p, length.out = n)
  lnmean = rep(lnmean, length.out = n)
  lnsd = rep(lnsd, length.out = n)
  u = rep(u, length.out = n)
  xi = rep(xi, length.out = n)
  
  pu = plnorm(u, lnmean, lnsd)
  if (is.logical(phiu)) {
    phiu = 1 - pu
  } else {
    phiu = rep(phiu, length.out = n)
  }
  phib = (1 - phiu) / pu
    
  sigmau = phiu / (phib * dlnorm(u, lnmean, lnsd))
  
  check.posparam(sigmau, allowvec = TRUE)
    
  qlognormgpd(p, lnmean, lnsd, u, sigmau, xi, phiu, lower.tail)
}

#' @export
#' @aliases lognormgpdcon dlognormgpdcon plognormgpdcon qlognormgpdcon rlognormgpdcon
#' @rdname  lognormgpdcon

# random number generation for log-normal bulk with GPD for upper tail
# with continuity at threshold
rlognormgpdcon <- function(n = 1, lnmean = 0, lnsd = 1, u = qlnorm(0.9, lnmean, lnsd),
  xi = 0, phiu = TRUE) {

  # Check properties of inputs
  check.n(n)
  check.param(lnmean, allowvec = TRUE)
  check.posparam(lnsd, allowvec = TRUE)
  check.posparam(u, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.phiu(phiu, allowvec = TRUE)

  n = check.inputn(c(n, length(lnmean), length(lnsd), length(u), length(xi), length(phiu)),
                   allowscalar = TRUE)

  if (any(xi == 1)) stop("shape cannot be 1")
  
  qlognormgpdcon(runif(n), lnmean, lnsd, u, xi, phiu)
}
