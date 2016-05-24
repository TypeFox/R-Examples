#' @name gngcon
#' 
#' @title Normal Bulk with GPD Upper and Lower Tails Extreme Value Mixture Model
#'  with Single Continuity Constraint at Thresholds
#'
#' @description Density, cumulative distribution function, quantile function and
#'   random number generation for the extreme value mixture model with normal
#'   for bulk distribution between the upper and lower thresholds with
#'   conditional GPD's for the two tails with continuity at the lower and upper thresholds.
#'   The parameters are the normal mean
#'   \code{nmean} and standard deviation \code{nsd}, lower tail (threshold \code{ul}, 
#'   GPD shape \code{xil} and tail fraction \code{phiul})
#'   and upper tail (threshold \code{ur}, GPD shape 
#'   \code{xiR} and tail fraction \code{phiuR}).
#'
#' @inheritParams gng
#' 
#' @details Extreme value mixture model combining normal distribution for the bulk
#' between the lower and upper thresholds and GPD for upper and lower tails with Continuity Constraints at the lower and upper threshold. The
#' user can pre-specify \code{phiul} and \code{phiur} permitting a parameterised
#' value for the lower and upper tail fraction respectively. Alternatively, when
#' \code{phiul=TRUE} or \code{phiur=TRUE} the corresponding tail fraction is
#' estimated as from the normal bulk model.
#' 
#' Notice that the tail fraction cannot be 0 or 1, and the sum of upper and lower tail
#' fractions \code{phiul+phiur<1}, so the lower threshold must be less than the upper, 
#' \code{ul<ur}.
#' 
#' The cumulative distribution function now has three components. The lower tail with 
#' tail fraction \eqn{\phi_{ul}} defined by the normal bulk model (\code{phiul=TRUE})
#' upto the lower threshold \eqn{x < u_l}:
#'   \deqn{F(x) = H(u_l) G_l(x).}
#' where \eqn{H(x)} is the normal cumulative distribution function (i.e. 
#' \code{pnorm(ur, nmean, nsd)}). The 
#' \eqn{G_l(X)} is the conditional GPD cumulative distribution function with negated
#' data and threshold, i.e. \code{dgpd(-x, -ul, sigmaul, xil, phiul)}. The normal
#' bulk model between the thresholds \eqn{u_l \le x \le u_r} given by:
#' \deqn{F(x) = H(x).}
#' Above the threshold \eqn{x > u_r} the usual conditional GPD:
#' \deqn{F(x) = H(u_r) + [1 - H(u_r)] G(x)}
#' where \eqn{G(X)}.
#' 
#' The cumulative distribution function for the pre-specified tail fractions 
#' \eqn{\phi_{ul}} and \eqn{\phi_{ur}} is more complicated.  The unconditional GPD
#' is used for the lower tail \eqn{x < u_l}:
#'   \deqn{F(x) = \phi_{ul} G_l(x).}
#' The normal bulk model between the thresholds \eqn{u_l \le x \le u_r} given by:
#' \deqn{F(x) = \phi_{ul}+ (1-\phi_{ul}-\phi_{ur}) (H(x) - H(u_l)) / (H(u_r) - H(u_l)).}
#' Above the threshold \eqn{x > u_r} the usual conditional GPD:
#' \deqn{F(x) = (1-\phi_{ur}) + \phi_{ur} G(x)}
#' Notice that these definitions are equivalent when \eqn{\phi_{ul} = H(u_l)} and
#' \eqn{\phi_{ur} = 1 - H(u_r)}.
#'  
#' The continuity constraint at \code{ur} means that:
#'  \deqn{\phi_{ur} g_r(x) = (1-\phi_{ul}-\phi_{ur}) h(u_l)/ (H(u_r) - H(u_l)).}
#' By rearrangement, the GPD scale parameter \code{sigmaur} is then:
#' \deqn{\sigma_ur = \phi_{ur} (H(u_r) - H(u_l))/ h(u_l) (1-\phi_{ul}-\phi_{ur}).}
#' where \eqn{h(x)}, \eqn{g_l(x)} and \eqn{g_r(x)} are the normal and conditional GPD
#' density functions for lower and upper tail respectively. 
#' In the special case of where the tail fraction is defined by the bulk model this reduces to
#' \deqn{\sigma_ur = [1-H(u_r)] / h(u_r)}.
#' 
#' The continuity constraint at \code{ul} means that:
#'  \deqn{\phi_{ul} g_l(x) = (1-\phi_{ul}-\phi_{ur}) h(u_l)/ (H(u_r) - H(u_l)).}
#' The GPD scale parameter \code{sigmaul} is replaced by:
#' \deqn{\sigma_ul = \phi_{ul} (H(u_r) - H(u_l))/ h(u_l) (1-\phi_{ul}-\phi_{ur}).}
#' In the special case of where the tail fraction is defined by the bulk model this reduces to
#' \deqn{\sigma_ul = H(u_l)/ h(u_l)}. 
#' 
#' See \code{\link[evmix:gpd]{gpd}} for details of GPD upper tail component, 
#' \code{\link[stats:Normal]{dnorm}} for details of normal bulk component,
#' \code{\link[evmix:normgpd]{dnormgpd}} for normal with GPD extreme value
#' mixture model and \code{\link[evmix:gng]{dgng}} for normal bulk with GPD 
#'  upper and lower tails extreme value mixture model.
#' 
#' @return \code{\link[evmix:gngcon]{dgngcon}} gives the density, 
#' \code{\link[evmix:gngcon]{pgngcon}} gives the cumulative distribution function,
#' \code{\link[evmix:gngcon]{qgngcon}} gives the quantile function and 
#' \code{\link[evmix:gngcon]{rgngcon}} gives a random sample.
#' 
#' @note All inputs are vectorised except \code{log} and \code{lower.tail}.
#' The main inputs (\code{x}, \code{p} or \code{q}) and parameters must be either
#' a scalar or a vector. If vectors are provided they must all be of the same length,
#' and the function will be evaluated for each element of vector. In the case of 
#' \code{\link[evmix:gngcon]{rgngcon}} any input vector must be of length \code{n}.
#' 
#' Default values are provided for all inputs, except for the fundamentals 
#' \code{x}, \code{q} and \code{p}. The default sample size for 
#' \code{\link[evmix:gngcon]{rgngcon}} is 1.
#' 
#' Missing (\code{NA}) and Not-a-Number (\code{NaN}) values in \code{x},
#' \code{p} and \code{q} are passed through as is and infinite values are set to
#' \code{NA}. None of these are not permitted for the parameters.
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
#' Zhao, X., Scarrott, C.J. Reale, M. and Oxley, L. (2010). Extreme value modelling
#' for forecasting the market crisis. Applied Financial Econometrics 20(1), 63-72.
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @seealso \code{\link[evmix:gpd]{gpd}} and \code{\link[stats:Normal]{dnorm}}
#' @aliases gngcon dgngcon pgngcon qgngcon rgngcon
#' @family  normgpd normgpdcon gng gngcon fnormgpd fnormgpdcon fgng fgngcon
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' par(mfrow = c(2, 2))
#' 
#' x = rgngcon(1000, phiul = 0.15, phiur = 0.15)
#' xx = seq(-6, 6, 0.01)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-6, 6))
#' lines(xx, dgngcon(xx, phiul = 0.15, phiur = 0.15))
#' 
#' # three tail behaviours
#' plot(xx, pgngcon(xx), type = "l")
#' lines(xx, pgngcon(xx, xil = 0.3, xir = 0.3), col = "red")
#' lines(xx, pgngcon(xx, xil = -0.3, xir = -0.3), col = "blue")
#' legend("topleft", paste("Symmetric xil=xir=",c(0, 0.3, -0.3)),
#'   col=c("black", "red", "blue"), lty = 1)
#' 
#' x = rgngcon(1000, xil = -0.3, phiul = 0.2, xir = 0.3, phiur = 0.2)
#' xx = seq(-6, 6, 0.01)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-6, 6))
#' lines(xx, dgngcon(xx, xil = -0.3, phiul = 0.2, xir = 0.3, phiur = 0.2))
#' 
#' plot(xx, dgngcon(xx, xil = -0.3, phiul = 0.2, xir = 0.3, phiur = 0.2), type = "l", ylim = c(0, 0.4))
#' lines(xx, dgngcon(xx, xil = -0.3, phiul = 0.3, xir = 0.3, phiur = 0.3), col = "red")
#' lines(xx, dgngcon(xx, xil = -0.3, phiul = TRUE, xir = 0.3, phiur = TRUE), col = "blue")
#' legend("topleft", c("phiul = phiur = 0.2", "phiul = phiur = 0.3", "Bulk Tail Fraction"),
#'   col=c("black", "red", "blue"), lty = 1)
#' }
#' 
NULL

#' @export
#' @aliases gngcon dgngcon pgngcon qgngcon rgngcon
#' @rdname  gngcon

# probability density function for normal bulk with GPD's for upper and lower tails
# with continuity at both thresholds
dgngcon <- function(x, nmean = 0, nsd = 1,
  ul = qnorm(0.1, nmean, nsd), xil = 0, phiul = TRUE, 
  ur = qnorm(0.9, nmean, nsd), xir = 0, phiur = TRUE, log = FALSE) {
  
  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(nmean, allowvec = TRUE)
  check.posparam(nsd, allowvec = TRUE)
  check.param(ul, allowvec = TRUE)
  check.param(xil, allowvec = TRUE)
  check.phiu(phiul, allowvec = TRUE)
  check.param(ur, allowvec = TRUE)
  check.param(xir, allowvec = TRUE)
  check.phiu(phiur, allowvec = TRUE)
  check.logic(log)

  n = check.inputn(c(length(x), length(nmean), length(nsd),
    length(ul), length(xil), length(phiul),
    length(ur), length(xir), length(phiur)), allowscalar = TRUE)

  if (any(is.infinite(x))) warning("infinite quantiles set to NA")

  x[is.infinite(x)] = NA # user will have to deal with infinite cases

  if (any(ul >= ur)) stop("lower threshold must be below upper threshold")
  
  if (!is.logical(phiul) & !is.logical(phiur)) {
    if (any((phiul + phiur) > 1)) stop("phiu + phiur must be less than 1")
  }
  
  x = rep(x, length.out = n)
  nmean = rep(nmean, length.out = n)
  nsd = rep(nsd, length.out = n)
  ul = rep(ul, length.out = n)
  xil = rep(xil, length.out = n)
  ur = rep(ur, length.out = n)
  xir = rep(xir, length.out = n)
  
  if (is.logical(phiul)) {
    phiul = pnorm(ul, nmean, nsd)
  } else {
    phiul = rep(phiul, length.out = n)
  }
  if (is.logical(phiur)) {
    phiur = 1 - pnorm(ur, nmean, nsd)
  } else {
    phiur = rep(phiur, length.out = n)
  }
  phib = (1 - phiul - phiur) / (pnorm(ur, nmean, nsd) - pnorm(ul, nmean, nsd))
  
  sigmaul = phiul / (phib * dnorm(ul, nmean, nsd))
  sigmaur = phiur / (phib * dnorm(ur, nmean, nsd)) 
  
  check.posparam(sigmaul, allowvec = TRUE)
  check.posparam(sigmaur, allowvec = TRUE)
    
  dgng(x, nmean, nsd, ul, sigmaul, xil, phiul, ur, sigmaur, xir, phiur, log)
}

#' @export
#' @aliases gngcon dgngcon pgngcon qgngcon rgngcon
#' @rdname  gngcon

# cumulative distribution function for normal bulk with GPD's for upper and lower tails
# with continuity at both thresholds
pgngcon <- function(q, nmean = 0, nsd = 1,
  ul = qnorm(0.1, nmean, nsd), xil = 0, phiul = TRUE, 
  ur = qnorm(0.9, nmean, nsd), xir = 0, phiur = TRUE, lower.tail = TRUE) {

  # Check properties of inputs
  check.quant(q, allowna = TRUE, allowinf = TRUE)
  check.param(nmean, allowvec = TRUE)
  check.posparam(nsd, allowvec = TRUE)
  check.param(ul, allowvec = TRUE)
  check.param(xil, allowvec = TRUE)
  check.phiu(phiul, allowvec = TRUE)
  check.param(ur, allowvec = TRUE)
  check.param(xir, allowvec = TRUE)
  check.phiu(phiur, allowvec = TRUE)
  check.logic(lower.tail)

  n = check.inputn(c(length(q), length(nmean), length(nsd),
    length(ul), length(xil), length(phiul),
    length(ur), length(xir), length(phiur)), allowscalar = TRUE)

  if (any(is.infinite(q))) warning("infinite quantiles set to NA")

  q[is.infinite(q)] = NA # user will have to deal with infinite cases

  if (any(ul >= ur)) stop("lower threshold must be below upper threshold")

  if (!is.logical(phiul) & !is.logical(phiur)) {
    if (any((phiul + phiur) > 1)) stop("phiu + phiur must be less than 1")
  }

  q = rep(q, length.out = n)
  nmean = rep(nmean, length.out = n)
  nsd = rep(nsd, length.out = n)
  ul = rep(ul, length.out = n)
  xil = rep(xil, length.out = n)
  ur = rep(ur, length.out = n)
  xir = rep(xir, length.out = n)
  
  if (is.logical(phiul)) {
    phiul = pnorm(ul, nmean, nsd)
  } else {
    phiul = rep(phiul, length.out = n)
  }
  if (is.logical(phiur)) {
    phiur = 1 - pnorm(ur, nmean, nsd)
  } else {
    phiur = rep(phiur, length.out = n)
  }
  phib = (1 - phiul - phiur) / (pnorm(ur, nmean, nsd) - pnorm(ul, nmean, nsd))
  
  sigmaul = phiul / (phib * dnorm(ul, nmean, nsd))
  sigmaur = phiur / (phib * dnorm(ur, nmean, nsd)) 
  
  check.posparam(sigmaul, allowvec = TRUE)
  check.posparam(sigmaur, allowvec = TRUE)
    
  pgng(q, nmean, nsd, ul, sigmaul, xil, phiul, ur, sigmaur, xir, phiur, lower.tail)
  
}

#' @export
#' @aliases gngcon dgngcon pgngcon qgngcon rgngcon
#' @rdname  gngcon

# inverse cumulative distribution function for normal bulk with GPD's for upper and lower tails
# with continuity at both thresholds
qgngcon <- function(p, nmean = 0, nsd = 1,
  ul = qnorm(0.1, nmean, nsd), xil = 0, phiul = TRUE, 
  ur = qnorm(0.9, nmean, nsd), xir = 0, phiur = TRUE, lower.tail = TRUE) {
  
  # Check properties of inputs
  check.prob(p, allowna = TRUE)
  check.param(nmean, allowvec = TRUE)
  check.posparam(nsd, allowvec = TRUE)
  check.param(ul, allowvec = TRUE)
  check.param(xil, allowvec = TRUE)
  check.phiu(phiul, allowvec = TRUE)
  check.param(ur, allowvec = TRUE)
  check.param(xir, allowvec = TRUE)
  check.phiu(phiur, allowvec = TRUE)
  check.logic(lower.tail)

  n = check.inputn(c(length(p), length(nmean), length(nsd),
    length(ul), length(xil), length(phiul),
    length(ur), length(xir), length(phiur)), allowscalar = TRUE)

  if (any(ul >= ur)) stop("lower threshold must be below upper threshold")

  if (!is.logical(phiul) & !is.logical(phiur)) {
    if (any((phiul + phiur) > 1)) stop("phiu + phiur must be less than 1")
  }
    
  p = rep(p, length.out = n)
  nmean = rep(nmean, length.out = n)
  nsd = rep(nsd, length.out = n)
  ul = rep(ul, length.out = n)
  xil = rep(xil, length.out = n)
  ur = rep(ur, length.out = n)
  xir = rep(xir, length.out = n)
  
  if (is.logical(phiul)) {
    phiul = pnorm(ul, nmean, nsd)
  } else {
    phiul = rep(phiul, length.out = n)
  }
  if (is.logical(phiur)) {
    phiur = 1 - pnorm(ur, nmean, nsd)
  } else {
    phiur = rep(phiur, length.out = n)
  }
  phib = (1 - phiul - phiur) / (pnorm(ur, nmean, nsd) - pnorm(ul, nmean, nsd))
  
  sigmaul = phiul / (phib * dnorm(ul, nmean, nsd))
  sigmaur = phiur / (phib * dnorm(ur, nmean, nsd)) 
  
  check.posparam(sigmaul, allowvec = TRUE)
  check.posparam(sigmaur, allowvec = TRUE)
    
  qgng(p, nmean, nsd, ul, sigmaul, xil, phiul, ur, sigmaur, xir, phiur, lower.tail)
}

#' @export
#' @aliases gngcon dgngcon pgngcon qgngcon rgngcon
#' @rdname  gngcon

# random number generation for normal bulk with GPD's for upper and lower tails
# with continuity at both thresholds
rgngcon <- function(n = 1, nmean = 0, nsd = 1,
  ul = qnorm(0.1, nmean, nsd), xil = 0, phiul = TRUE, 
  ur = qnorm(0.9, nmean, nsd), xir = 0, phiur = TRUE) {
  
  # Check properties of inputs
  check.n(n)
  check.param(nmean, allowvec = TRUE)
  check.posparam(nsd, allowvec = TRUE)
  check.param(ul, allowvec = TRUE)
  check.param(xil, allowvec = TRUE)
  check.phiu(phiul, allowvec = TRUE)
  check.param(ur, allowvec = TRUE)
  check.param(xir, allowvec = TRUE)
  check.phiu(phiur, allowvec = TRUE)

  n = check.inputn(c(n, length(nmean), length(nsd),
    length(ul), length(xil), length(phiul),
    length(ur), length(xir), length(phiur)), allowscalar = TRUE)

  if (any(xil == 1) | any(xir == 1)) stop("shape cannot be 1")
  
  qgngcon(runif(n), nmean, nsd, ul, xil, phiul, ur, xir, phiur)
}
