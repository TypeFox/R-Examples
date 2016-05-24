#' @name psden
#' 
#' @title P-Splines probability density function
#'
#' @description Density, cumulative distribution function, quantile function and
#'   random number generation for the P-splines density estimate. B-spline coefficients
#'   can be result from Poisson regression with log or identity link.
#'
#' @inheritParams gpd
#' @param beta          vector of B-spline coefficients (required)
#' @param nbinwidth     scaling to convert count frequency into proper density
#' @param xrange        vector of minimum and maximum of B-spline (support of density)
#' @param nseg          number of segments between knots
#' @param degree        degree of B-splines (0 is constant, 1 is linear, etc.)
#' @param design.knots  spline knots for splineDesign function
#' 
#' @details P-spline density estimate using B-splines with given coefficients. B-splines
#' knots can be specified using \code{design.knots} or regularly spaced knots can be specified
#' using \code{xrange}, \code{nseg} and \code{deg}. No default knots are provided.
#' 
#' If regularly spaced knots are specified using \code{xrange}, \code{nseg} and \code{deg},
#' then B-splines which are shifted/spliced versions of each other are defined (i.e. not natural B-splines)
#' which is consistent with definition of Eilers and Marx, the masters of P-splines.
#' 
#' The \code{\link[splines:splineDesign]{splineDesign}} function is used to calculate the B-splines, which 
#' intakes knot locations as \code{design.knots}. As such the \code{design.knots} are not the knots in
#' their usual sense (e.g. to cover [0, 100] with 10 segments the usual knots would be \eqn{0, 10, \ldots, 100}).
#' The \code{design.knots} must be extended by the \code{degree}, so for \code{degree = 2} the
#' \code{design.knots = seq(-20, 120, 10)}.
#' 
#' Further, if the user wants natural B-splines then these can be specified using the
#' \code{design.knots}, with replicated knots at each bounday according to the degree. To continue the 
#' above example, for \code{degree = 2} the \code{design.knots = c(rep(0, 2), seq(0, 100, 10), rep(100, 2))}. 
#' 
#' If both the \code{design.knots} and other knot specification are provided, then the former are
#' used by default. Default values for only the \code{degree} and \code{nseg} are provided, all the other
#' P-spline inputs must be provided. Notice that the \code{order} and \code{lambda} penalty are not needed
#' as these are encapsulated in the inference for the B-spline coefficients.
#' 
#' Poisson regression is typically used for estimating the B-spline coefficients, using maximum likelihood
#' estimation (via iterative re-weighted least squares). A log-link function is usually used and as such the 
#' \code{beta} coefficients are on a log-scale, and the density needs to be exponentiated. However, an
#' identity link may be (carefully) used and then these coefficients are on the usual scale.
#' 
#' The \code{beta} coefficients are estimated using a particular sample (size) and histogram bin-width, using 
#' Poisson regression. Thus to
#' convert the predicted counts into a proper density it needs to be rescaled by dividing by \eqn{n * binwidth}.
#' If \code{nbinwidth=NULL} is not provided then a crude approximate scaling is used by normalising the density
#' to be proper. The renormalisation requires numerical integration, which is
#' computationally intensive and so best avoided wherever possible.
#' 
#' Checks of the consistency of the \code{xrange}, \code{degree} and \code{nseg} and \code{design.knots} are made,
#' with the values implied by the \code{design.knots} used by default to replace any incorrect values. These
#' replacements are made for notational efficiency for users.
#' 
#' An inversion sampler is used for random number generation which also rather
#' inefficient, as it could be carried out more efficiently using a mixture representation.
#' 
#' The quantile function is rather complicated as there is no closed form solution,
#' so is obtained by numerical approximation of the inverse cumulative distribution function
#' \eqn{P(X \le q) = p} to find \eqn{q}. The quantile function 
#' \code{\link[evmix:psden]{qpsden}} evaluates the P-splines cumulative distribution
#' function over the \code{xrange}. A sequence of values
#' of length fifty times the number of knots (with a minimum of 1000) is first
#' calculated. Spline based interpolation using \code{\link[stats:splinefun]{splinefun}},
#' with default \code{monoH.FC} method, is then used to approximate the quantile
#' function. This is a similar approach to that taken
#' by Matt Wand in the \code{\link[ks:kde.1d]{qkde}} in the \code{\link[ks:kde.1d]{ks}} package.
#' 
#' @return \code{\link[evmix:psden]{dpsden}} gives the density, 
#' \code{\link[evmix:psden]{ppsden}} gives the cumulative distribution function,
#' \code{\link[evmix:psden]{qpsden}} gives the quantile function and 
#' \code{\link[evmix:psden]{rpsden}} gives a random sample.
#' 
#' @note Unlike most of the other extreme value mixture model functions the 
#'   \code{\link[evmix:psden]{psden}} functions have not been vectorised as
#'   this is not appropriate. The main inputs (\code{x}, \code{p} or \code{q})
#'   must be either a scalar or a vector, which also define the output length.
#' 
#' Default values are provided for P-spline inputs of \code{degree} and \code{nseg} only, 
#' but all others must be provided by the user.
#' The default sample size for \code{\link[evmix:psden]{rpsden}} is 1.
#' 
#' Missing (\code{NA}) and Not-a-Number (\code{NaN}) values in \code{x},
#' \code{p} and \code{q} are passed through as is and infinite values are set to
#' \code{NA}. None of these are not permitted for the parameters.
#' 
#' Error checking of the inputs (e.g. invalid probabilities) is carried out and
#' will either stop or give warning message as appropriate.
#' 
#' @references
#' \url{http://en.wikipedia.org/wiki/B-spline}
#' 
#' \url{http://www.stat.lsu.edu/faculty/marx/}
#' 
#' Eilers, P.H.C. and Marx, B.D. (1996). Flexible smoothing with B-splines and penalties.
#' Statistical Science 11(2), 89-121.
#' 
#' @author Alfadino Akbar and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}.
#'
#' @seealso \code{\link[splines:splineDesign]{splineDesign}}.
#' 
#' @aliases psden dpsden ppsden qpsden rpsden
#' @family  psden fpsden
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
#' fit = fpsden(x, lambdaseq = 10^seq(-5, 5, 0.25), breaks = breaks,
#'              xrange = c(-4, 4), nseg = 10, degree = 3, ord = 2)
#' psdensity = exp(fit$bsplines %*% fit$mle)
#' 
#' hist(x, freq = FALSE, breaks = seq(-4, 4, length.out=101), xlim = c(-6, 6))
#' lines(xx, y, col = "black") # true density
#' 
#' # P-splines density from dpsden function
#' with(fit, lines(xx, dpsden(xx, beta, nbinwidth, design = design.knots), lwd = 2, col = "blue"))
#'
#' legend("topright", c("True Density","P-spline density"), col=c("black", "blue"), lty = 1)
#' 
#' # plot B-splines
#' par(mfrow = c(2, 1))
#' with(fit, matplot(mids, as.matrix(bsplines), type = "l", lty = 1))
#' 
#' # Natural B-splines
#' knots = with(fit, seq(xrange[1], xrange[2], length.out = nseg + 1))
#' natural.knots = with(fit, c(rep(xrange[1], degree), knots, rep(xrange[2], degree)))
#' naturalb = splineDesign(natural.knots, fit$mids, ord = fit$degree + 1, outer.ok = TRUE)
#' with(fit, matplot(mids, naturalb, type = "l", lty = 1))
#'
#' # Compare knot specifications
#' rbind(fit$design.knots, natural.knots)
#' 
#' # User can use natural B-splines if design.knots are specified manually
#' natural.fit = fpsden(x, lambdaseq = 10^seq(-5, 5, 0.25), breaks = breaks,
#'              design.knots = natural.knots, nseg = 10, degree = 3, ord = 2)
#' psdensity = with(natural.fit, exp(bsplines %*% mle))
#' 
#' par(mfrow = c(1, 1))
#' hist(x, freq = FALSE, breaks = seq(-4, 4, length.out=101), xlim = c(-6, 6))
#' lines(xx, y, col = "black") # true density
#' 
#' # check density against dpsden function
#' with(fit, lines(xx, dpsden(xx, beta, nbinwidth, design = design.knots), lwd = 2, col = "blue"))
#' with(natural.fit, lines(xx, dpsden(xx, beta, nbinwidth, design = design.knots),
#'                         lwd = 2, col = "red", lty = 2))
#'
#' legend("topright", c("True Density", "Eilers and Marx B-splines", "Natural B-splines"),
#'    col=c("black", "blue", "red"), lty = c(1, 1, 2))
#' }
#' 
NULL

#' @export
#' @aliases psden dpsden ppsden qpsden rpsden
#' @rdname  psden

# density function for P-splines density estimator
dpsden <- function(x, beta = NULL, nbinwidth = NULL, xrange = NULL, nseg = 10, degree = 3, design.knots = NULL,
                   log = FALSE) {

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(beta, allowvec = TRUE)
  check.posparam(nbinwidth, allownull = TRUE)
  check.param(xrange, allowvec = TRUE, allownull = TRUE)
  check.n(nseg)
  check.n(degree, allowzero = TRUE)
  check.param(design.knots, allowvec = TRUE, allownull = TRUE)
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

  d = x # pass through NA/NaN as entered

  whichok = which(is.finite(x))

  # Only interpolate with B-splines, no extrapolation
  xint = x[whichok]
  d[whichok] = -Inf
  whichint = which((xint >= xrange[1]) & (xint <= xrange[2]))
    
  bsplines = splineDesign(design.knots, xint[whichint], ord = degree + 1)

  d[whichok[whichint]] = bsplines %*% beta - log(nbinwidth) # equivalent to dividing by nbinwidth on original scale
    
  if (!log) d = exp(d)
  
  return(d)
}

#' @export
#' @aliases psden dpsden ppsden qpsden rpsden
#' @rdname  psden

# cumulative distribution function for P-splines density estimator
ppsden <- function(q, beta = NULL, nbinwidth = NULL, xrange = NULL, nseg = 10, degree = 3, design.knots = NULL,
                   lower.tail = TRUE) {

  # Check properties of inputs
  check.quant(q, allowna = TRUE, allowinf = TRUE)
  check.param(beta, allowvec = TRUE)
  check.posparam(nbinwidth, allownull = TRUE)
  check.param(xrange, allowvec = TRUE, allownull = TRUE) 
  check.n(nseg)
  check.n(degree, allowzero = TRUE)
  check.param(design.knots, allowvec = TRUE, allownull = TRUE) 
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

  p = q # pass through NA/NaN as entered
  
  whichok = which(is.finite(q))
  
  # Only interpolate with B-splines, no extrapolation
  qint = q[whichok]
  p[whichok] = ifelse(qint < xrange[1], 0, ifelse(qint > xrange[2], 1, NA))
  whichint = which((qint >= xrange[1]) & (qint <= xrange[2]))

  # numerical integration required
  intpsden <- function(x, beta, design.knots, degree) {
    psdenint = try(integrate(pscounts, lower = xrange[1], upper = x,
      beta = beta, design.knots = design.knots, degree = degree,
      subdivisions = 10000, rel.tol = 1e-9, stop.on.error = FALSE))
    
    if (inherits(psdenint, "try-error")) {
      psdenint$value = NA
      warning("failed to numerically evaluate cdf of P-spline density estimate")
    }
    
    psdenint$value
  }
  
  p[whichok[whichint]] = sapply(qint[whichint], intpsden,
                                beta = beta, design.knots = design.knots, degree = degree) / nbinwidth

  # sometimes due to numerical errors p>1 or p<0
  p = pmax(pmin(p, 1), 0)
  
  if (!lower.tail) p = 1 - p

  p
}

#' @export
#' @aliases psden dpsden ppsden qpsden rpsden
#' @rdname  psden

# inverse cumulative distribution function for P-splines density estimator
qpsden <- function(p, beta = NULL, nbinwidth = NULL, xrange = NULL, nseg = 10, degree = 3, design.knots = NULL,
                   lower.tail = TRUE) {

  # Check properties of inputs
  check.prob(p, allowna = TRUE)
  check.param(beta, allowvec = TRUE)
  check.posparam(nbinwidth, allownull = TRUE)
  check.param(xrange, allowvec = TRUE, allownull = TRUE) 
  check.n(nseg)
  check.n(degree, allowzero = TRUE)
  check.param(design.knots, allowvec = TRUE, allownull = TRUE) 
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

  if (!lower.tail) p = 1 - p

  q = p
  
  # obtain quantile function but interpolation between P-spline based CDF estimates
  # - CDF is reasonably quick to estimate, interpolation algorithms also quite fast
  # - an alternative solution is to numerical solve to find quantile, but is much slower
  
  qk = seq(xrange[1], xrange[2], length.out = min(50 * length(design.knots), 1000))

  pk = ppsden(qk, beta = beta, nbinwidth = nbinwidth,
              xrange = xrange, nseg = nseg, degree = degree, design.knots = design.knots)

  qfun = splinefun(x = pk, y = qk) # spline based inverse function

  whichok = which(is.finite(q))
  q[whichok] = qfun(p[whichok])

  q
}


#' @export
#' @aliases psden dpsden ppsden qpsden rpsden
#' @rdname  psden

# cumulative distribution function for P-splines density estimator
rpsden <- function(n = 1, beta = NULL, nbinwidth = NULL, xrange = NULL, nseg = 10, degree = 3,
                   design.knots = NULL) {

  # Check properties of inputs
  check.n(n)
  check.param(beta, allowvec = TRUE)
  check.posparam(nbinwidth, allownull = TRUE)
  check.param(xrange, allowvec = TRUE, allownull = TRUE) 
  check.n(nseg)
  check.n(degree, allowzero = TRUE)
  check.param(design.knots, allowvec = TRUE, allownull = TRUE) 

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

  # Inversion sampling (inefficient, but works)
  qpsden(runif(n), beta, nbinwidth, xrange, nseg, degree, design.knots)
}
