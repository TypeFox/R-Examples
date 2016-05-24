#' @name psdengpd
#' 
#' @title P-Splines Density Estimate and GPD Tail Extreme Value Mixture Model
#'
#' @description Density, cumulative distribution function, quantile function and
#'   random number generation for the extreme value mixture model with P-splines density estimate for bulk
#'   distribution upto the threshold and conditional GPD above threshold. The parameters
#'   are the B-spline coefficients \code{beta} (and associated features), threshold \code{u}
#'   GPD scale \code{sigmau} and shape \code{xi} and tail fraction \code{phiu}.
#'
#' @inheritParams normgpd
#' @inheritParams psden
#' @inheritParams gpd
#' 
#' @details Extreme value mixture model combining P-splines density estimate for the bulk
#' below the threshold and GPD for upper tail.
#' 
#' The user can pre-specify \code{phiu} 
#' permitting a parameterised value for the tail fraction \eqn{\phi_u}. Alternatively, when
#' \code{phiu=TRUE} the tail fraction is estimated as the tail fraction from the
#' KDE bulk model.
#' 
#' The cumulative distribution function with tail fraction \eqn{\phi_u} defined by the
#' upper tail fraction of the P-splines density estimate (\code{phiu=TRUE}), upto the 
#' threshold \eqn{x \le u}, given by:
#' \deqn{F(x) = H(x)}
#' and above the threshold \eqn{x > u}:
#' \deqn{F(x) = H(u) + [1 - H(u)] G(x)}
#' where \eqn{H(x)} and \eqn{G(X)} are the P-splines density estimate and conditional GPD
#' cumulative distribution functions respectively.
#' 
#' The cumulative distribution function for pre-specified \eqn{\phi_u}, upto the
#' threshold \eqn{x \le u}, is given by:
#' \deqn{F(x) = (1 - \phi_u) H(x)/H(u)}
#' and above the threshold \eqn{x > u}:
#' \deqn{F(x) = \phi_u + [1 - \phi_u] G(x)}
#' Notice that these definitions are equivalent when \eqn{\phi_u = 1 - H(u)}.
#' 
#' See \code{\link[evmix:gpd]{gpd}} for details of GPD upper tail component. 
#' The specification of the underlying B-splines and the P-splines density estimator
#' are discussed in the \code{\link[evmix:psden]{psden}} function help.
#' 
#' @return \code{\link[evmix:psdengpd]{dpsdengpd}} gives the density, 
#' \code{\link[evmix:psdengpd]{ppsdengpd}} gives the cumulative distribution function,
#' \code{\link[evmix:psdengpd]{qpsdengpd}} gives the quantile function and 
#' \code{\link[evmix:psdengpd]{rpsdengpd}} gives a random sample.
#' 
#' @note Unlike most of the other extreme value mixture model functions the 
#' \code{\link[evmix:psdengpd]{psdengpd}} functions have not been vectorised as
#' this is not appropriate. The main inputs (\code{x}, \code{p} or \code{q})
#' must be either a scalar or a vector, which also define the output length.
#' The B-splines coefficients \code{beta} and knots \code{design.knots} are vectors.
#' 
#' Default values are provided for P-spline inputs of \code{degree} and \code{nseg} only, 
#' but all others must be provided by the user. The default sample size for
#' \code{\link[evmix:psdengpd]{rpsdengpd}} is 1.
#' 
#' Missing (\code{NA}) and Not-a-Number (\code{NaN}) values in \code{x},
#' \code{p} and \code{q} are passed through as is and infinite values are set to
#' \code{NA}. None of these are permitted for the parameters/B-spline criteria.
#' 
#' Due to symmetry, the lower tail can be described by GPD by negating the quantiles. 
#' 
#' Error checking of the inputs (e.g. invalid probabilities) is carried out and
#' will either stop or give warning message as appropriate.
#' 
#' @references
#' \url{http://en.wikipedia.org/wiki/B-spline}
#' 
#' \url{http://www.stat.lsu.edu/faculty/marx/}
#' 
#' \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
#' 
#' Scarrott, C.J. and MacDonald, A. (2012). A review of extreme value
#' threshold estimation and uncertainty quantification. REVSTAT - Statistical
#' Journal 10(1), 33-59. Available from \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}
#' 
#' Eilers, P.H.C. and Marx, B.D. (1996). Flexible smoothing with B-splines and penalties.
#' Statistical Science 11(2), 89-121.
#' 
#' @author Alfadino Akbar and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}.
#'
#' @seealso \code{\link[evmix:psden]{psden}} and \code{\link[evmix:fpsden]{fpsden}}.
#' 
#' @aliases psdengpd dpsdengpd ppsdengpd qpsdengpd rpsdengpd
#' @family  psden psdengpd fpsden fpsdengpd
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' par(mfrow = c(1, 1))
#' 
#' x = rnorm(1000)
#' xx = seq(-6, 6, 0.01)
#' y = dnorm(xx)
#' 
#' # Plenty of histogram bins (100)
#' breaks = seq(-4, 4, length.out=101)
#' 
#' # P-spline fitting with cubic B-splines, 2nd order penalty and 8 internal segments
#' # CV search for penalty coefficient. 
#' fit = fpsdengpd(x, lambdaseq = 10^seq(-5, 5, 0.25), breaks = breaks,
#'              xrange = c(-4, 4), nseg = 10, degree = 3, ord = 2)
#' hist(x, freq = FALSE, breaks = seq(-4, 4, length.out=101), xlim = c(-6, 6))
#' 
#' # P-splines only
#' with(fit, lines(xx, dpsden(xx, beta, nbinwidth, design = design.knots), lwd = 2, col = "blue"))
#'
#' # P-splines+GPD
#' with(fit, lines(xx, dpsdengpd(xx, beta, nbinwidth, design = design.knots, 
#'    u = u, sigmau = sigmau, xi = xi, phiu = phiu), lwd = 2, col = "red"))
#' abline(v = fit$u, col = "red")
#' 
#' legend("topleft", c("True Density","P-spline density", "P-spline+GPD"),
#'  col=c("black", "blue", "red"), lty = 1)
#' }
#' 
NULL

#' @export
#' @aliases psdengpd dpsdengpd ppsdengpd qpsdengpd rpsdengpd
#' @rdname  psdengpd

# probability density function for P-splines density estimate for the bulk
# distribution upto the threshold and conditional GPD above threshold
dpsdengpd <- function(x, beta = NULL, nbinwidth = NULL, xrange = NULL, nseg = 10, degree = 3,
                      u = NULL, sigmau = NULL, xi = 0, phiu = TRUE, design.knots = NULL, log = FALSE) {
  
  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(beta, allowvec = TRUE)
  check.posparam(nbinwidth, allownull = TRUE)
  check.param(xrange, allowvec = TRUE, allownull = TRUE)
  check.n(nseg)
  check.n(degree, allowzero = TRUE)
  check.param(design.knots, allowvec = TRUE, allownull = TRUE)
  check.param(u, allownull = TRUE)
  check.posparam(sigmau, allownull = TRUE)
  check.param(xi)
  check.phiu(phiu)
  check.logic(log)

  if (any(is.infinite(x))) warning("infinite quantiles set to NA")

  x[is.infinite(x)] = NA # user will have to deal with infinite cases

  checked.knots = check.design.knots(beta, xrange, nseg, degree, design.knots)
  xrange = checked.knots$xrange
  nseg = checked.knots$nseg
  degree = checked.knots$degree
  design.knots = checked.knots$design.knots

  # If constant for rescaling counts to density is not given, then density is renormalised to be proper
  if (is.null(nbinwidth)) {
    pscountint = try(integrate(pscounts, lower = xrange[1], upper = xrange[2],
      beta = beta, design.knots = design.knots, degree = degree,
      subdivisions = 10000, rel.tol = 1e-9, stop.on.error = FALSE))

    if (inherits(pscountint, "try-error"))
      stop("failed to numerically evaluate cdf of P-spline density estimate, provide nbinwidth input")
    nbinwidth = pscountint$value
  }
    
  # default for u and sigmau
  if (is.null(u)) {
    u = qpsden(0.9, beta, nbinwidth, xrange, nseg, degree, design.knots)
    check.param(u)
  }
  if (is.null(sigmau)) {
    sigmau = diff(qpsden(c(0.975, 0.85), beta, nbinwidth, xrange, nseg, degree, design.knots)) # if normal then approx sigma
    check.posparam(sigmau)
  }
  
  check.inputn(c(length(u), length(sigmau), length(xi), length(phiu)), allowscalar = TRUE) # scalar only
     
  pu = ppsden(u, beta, nbinwidth, xrange, nseg, degree, design.knots)
  if (is.logical(phiu)) {
    phiu = 1 - pu
  } else {
    phiu = phiu
  }
  phib = (1 - phiu) / pu

  d = x # pass through NA/NaN as entered

  whichb = which(x <= u)
  nb = length(whichb)
  whichu = which(x > u)
  nu = length(whichu)

  if (nb > 0) d[whichb] = log(phib) + dpsden(x[whichb], beta, nbinwidth,
                                             xrange, nseg, degree, design.knots, log = TRUE)

  if (nu > 0) d[whichu] = log(phiu) + dgpd(x[whichu], u, sigmau, xi, log = TRUE)
  
  if (!log) d = exp(d)
  
  d
}

#' @export
#' @aliases psdengpd dpsdengpd ppsdengpd qpsdengpd rpsdengpd
#' @rdname  psdengpd

# cumulative distribution function for P-splines density estimate for the bulk
# distribution upto the threshold and conditional GPD above threshold.
ppsdengpd <- function(q, beta = NULL, nbinwidth = NULL, xrange = NULL, nseg = 10, degree = 3,
                      u = NULL, sigmau = NULL, xi = 0, phiu = TRUE, design.knots = NULL, lower.tail = TRUE) {
  
  # Check properties of inputs
  check.quant(q, allowna = TRUE, allowinf = TRUE)
  check.param(beta, allowvec = TRUE)
  check.posparam(nbinwidth, allownull = TRUE)
  check.param(xrange, allowvec = TRUE, allownull = TRUE)
  check.n(nseg)
  check.n(degree, allowzero = TRUE)
  check.param(design.knots, allowvec = TRUE, allownull = TRUE)
  check.param(u, allownull = TRUE)
  check.posparam(sigmau, allownull = TRUE)
  check.param(xi)
  check.phiu(phiu)
  check.logic(lower.tail)

  if (any(is.infinite(q))) warning("infinite quantiles set to NA")

  q[is.infinite(q)] = NA # user will have to deal with infinite cases

  checked.knots = check.design.knots(beta, xrange, nseg, degree, design.knots)
  xrange = checked.knots$xrange
  nseg = checked.knots$nseg
  degree = checked.knots$degree
  design.knots = checked.knots$design.knots

  # If constant for rescaling counts to density is not given, then density is renormalised to be proper
  if (is.null(nbinwidth)) {
    pscountint = try(integrate(pscounts, lower = xrange[1], upper = xrange[2],
      beta = beta, design.knots = design.knots, degree = degree,
      subdivisions = 10000, rel.tol = 1e-9, stop.on.error = FALSE))

    if (inherits(pscountint, "try-error"))
      stop("failed to numerically evaluate cdf of P-spline density estimate, provide nbinwidth input")

    nbinwidth = pscountint$value
  }
    
  # default for u and sigmau
  if (is.null(u)) {
    u = qpsden(0.9, beta, nbinwidth, xrange, nseg, degree, design.knots)
    check.param(u)
  }
  if (is.null(sigmau)) {
    sigmau = diff(qpsden(c(0.975, 0.85), beta, nbinwidth, xrange, nseg, degree, design.knots)) # if normal then approx sigma
    check.posparam(sigmau)
  }

  check.inputn(c(length(u), length(sigmau), length(xi), length(phiu)), allowscalar = TRUE) # scalar only
  
  pu = ppsden(u, beta, nbinwidth, xrange, nseg, degree, design.knots)
  if (is.logical(phiu)) {
    phiu = 1 - pu
  } else {
    phiu = phiu
  }
  phib = (1 - phiu) / pu
    
  p = q # pass through NA/NaN as entered
  
  whichb = which(q <= u)
  nb = length(whichb)
  whichu = which(q > u)
  nu = length(whichu)
  
  if (nb > 0) p[whichb] = phib*ppsden(q[whichb], beta, nbinwidth, xrange, nseg, degree, design.knots)
  if (nu > 0) p[whichu] = (1 - phiu) + phiu*pgpd(q[whichu], u, sigmau, xi)
  
  if (!lower.tail) p = 1 - p
  
  p
}

#' @export
#' @aliases psdengpd dpsdengpd ppsdengpd qpsdengpd rpsdengpd
#' @rdname  psdengpd

# inverse cumulative distribution function for P-splines density estimate for the bulk
# distribution upto the threshold and conditional GPD above threshold.
qpsdengpd <- function(p, beta = NULL, nbinwidth = NULL, xrange = NULL, nseg = 10, degree = 3,
                      u = NULL, sigmau = NULL, xi = 0, phiu = TRUE, design.knots = NULL, lower.tail = TRUE) {
  
  # Check properties of inputs
  check.prob(p, allowna = TRUE)
  check.param(beta, allowvec = TRUE)
  check.posparam(nbinwidth, allownull = TRUE)
  check.param(xrange, allowvec = TRUE, allownull = TRUE)
  check.n(nseg)
  check.n(degree, allowzero = TRUE)
  check.param(design.knots, allowvec = TRUE, allownull = TRUE)
  check.param(u, allownull = TRUE)
  check.posparam(sigmau, allownull = TRUE)
  check.param(xi)
  check.phiu(phiu)
  check.logic(lower.tail)

  checked.knots = check.design.knots(beta, xrange, nseg, degree, design.knots)
  xrange = checked.knots$xrange
  nseg = checked.knots$nseg
  degree = checked.knots$degree
  design.knots = checked.knots$design.knots

  # If constant for rescaling counts to density is not given, then density is renormalised to be proper
  if (is.null(nbinwidth)) {
    pscountint = try(integrate(pscounts, lower = xrange[1], upper = xrange[2],
      beta = beta, design.knots = design.knots, degree = degree,
      subdivisions = 10000, rel.tol = 1e-9, stop.on.error = FALSE))

    if (inherits(pscountint, "try-error"))
      stop("failed to numerically evaluate cdf of P-spline density estimate, provide nbinwidth input")

    nbinwidth = pscountint$value
  }
    
  # default for u and sigmau
  if (is.null(u)) {
    u = qpsden(0.9, beta, nbinwidth, xrange, nseg, degree, design.knots)
    check.param(u)
  }
  if (is.null(sigmau)) {
    sigmau = diff(qpsden(c(0.975, 0.85), beta, nbinwidth, xrange, nseg, degree, design.knots)) # if normal then approx sigma
    check.posparam(sigmau)
  }

  check.inputn(c(length(u), length(sigmau), length(xi), length(phiu)), allowscalar = TRUE) # scalar only
  
  if (!lower.tail) p = 1 - p

  pu = ppsden(u, beta, nbinwidth, xrange, nseg, degree, design.knots)
  if (is.logical(phiu)) {
    phiu = 1 - pu
  } else {
    phiu = phiu
  }
  phib = (1 - phiu) / pu
  
  q = p # pass through NA/NaN as entered

  whichb = which(p <= (1 - phiu))
  nb = length(whichb)
  whichu = which(p > (1 - phiu))
  nu = length(whichu)
  
  if (nb > 0) q[whichb] = qpsden(p[whichb] / phib, beta, nbinwidth, xrange, nseg, degree, design.knots)
  if (nu > 0) q[whichu] = qgpd(p[whichu], u, sigmau, xi, phiu)
   
  q
}

#' @export
#' @aliases psdengpd dpsdengpd ppsdengpd qpsdengpd rpsdengpd
#' @rdname  psdengpd

# random number generation for P-splines density estimate for the bulk
# distribution upto the threshold and conditional GPD above threshold.
rpsdengpd <- function(n = 1, beta = NULL, nbinwidth = NULL, xrange = NULL, nseg = 10, degree = 3,
                      u = NULL, sigmau = NULL, xi = 0, phiu = TRUE, design.knots = NULL) {
  
  # Check properties of inputs
  check.n(n)
  check.param(beta, allowvec = TRUE)
  check.posparam(nbinwidth, allownull = TRUE)
  check.param(xrange, allowvec = TRUE, allownull = TRUE)
  check.n(nseg)
  check.n(degree, allowzero = TRUE)
  check.param(design.knots, allowvec = TRUE, allownull = TRUE)
  check.param(u, allownull = TRUE)
  check.posparam(sigmau, allownull = TRUE)
  check.param(xi)
  check.phiu(phiu)

  checked.knots = check.design.knots(beta, xrange, nseg, degree, design.knots)
  xrange = checked.knots$xrange
  nseg = checked.knots$nseg
  degree = checked.knots$degree
  design.knots = checked.knots$design.knots

  # If constant for rescaling counts to density is not given, then density is renormalised to be proper
  if (is.null(nbinwidth)) {
    pscountint = try(integrate(pscounts, lower = xrange[1], upper = xrange[2],
      beta = beta, design.knots = design.knots, degree = degree,
      subdivisions = 10000, rel.tol = 1e-9, stop.on.error = FALSE))

    if (inherits(pscountint, "try-error"))
      stop("failed to numerically evaluate cdf of P-spline density estimate, provide nbinwidth input")

    nbinwidth = pscountint$value
  }
    
  # default for u and sigmau
  if (is.null(u)) {
    u = qpsden(0.9, beta, nbinwidth, xrange, nseg, degree, design.knots)
    check.param(u)
  }
  if (is.null(sigmau)) {
    sigmau = diff(qpsden(c(0.975, 0.85), beta, nbinwidth, xrange, nseg, degree, design.knots)) # if normal then approx sigma
    check.posparam(sigmau)
  }

  check.inputn(c(length(u), length(sigmau), length(xi), length(phiu)), allowscalar = TRUE) # scalar only
    
  if (any(xi == 1)) stop("shape cannot be 1")
  
  qpsdengpd(runif(n), beta, nbinwidth, xrange, nseg, degree, u, sigmau, xi, phiu, design.knots)
}
