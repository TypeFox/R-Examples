#' @export
#' 
#' @title MLE Fitting of P-splines Density Estimator
#'
#' @description Maximum likelihood estimation for P-splines density estimation. Histogram binning
#' produces frequency counts, which are modelled by constrained B-splines in a Poisson regression. A penalty
#' based on differences in the sequences B-spline coefficients is used to smooth/interpolate the counts.
#' Iterated weighted least squares (IWLS) for a mixed model representation of the P-splines regression,
#' conditional on a particular penalty coefficient, is used for estimating the B-spline coefficients.
#' Leave-one-out cross-validation deviances are available for estimation of the penalty coefficient.
#' 
#' @param lambdaseq   vector of \eqn{\lambda}'s (or scalar) to be considered in profile likelihood. Required.
#' @param ord         order of difference used in the penalty term
#' @param lambda      penalty coefficient
#' @param breaks      histogram breaks (as in \code{\link[graphics:hist]{hist}} function)
#' @param counts      counts from histogram binning
#' @param bsplines    matrix of B-splines
#' @inheritParams     psden
#' @inheritParams     fgpd
#' 
#' @details The P-splines density estimator is fitted using maximum likelihood estimation, following
#' the approach of Eilers and Marx (1996). Histogram binning produces frequency counts, which are
#' modelled by constrained B-splines in a Poisson regression. A penalty
#' based on differences in the sequences B-spline coefficients is used to smooth/interpolate the counts.
#'
#' The B-splines are defined as in Eiler and Marx (1996), so that those are meet the boundary are simply
#' shifted and truncated version of the internal B-splines. No renormalisation is carried out. They are not
#' "natural" B-spline which are also commonly in use. Note that atural B-splines can be obtained by suitable
#' linear combinations of these B-splines. Hence, in practice there is little difference in the fit obtained
#' from either B-spline definition, even with the penalty constraining the coefficients. If the user desires
#' they can force the use of natural B-splines, by prior specification of the \code{design.knots}
#' with appropriate replication of the boundaries, see \code{\link[evmix:psden]{dpsden}}.
#' 
#' Iterated weighted least squares (IWLS) for a mixed model representation of the P-splines regression,
#' conditional on a particular penalty coefficient, is used for estimating the B-spline coefficients which
#' is equivalent to maximum likelihood estimation. Leave-one-out cross-validation deviances are available
#' for estimation of the penalty coefficient.
#' 
#' The parameter vector is the B-spline coefficients \code{beta}, no matter whether the penalty coefficient is
#' fixed or estimated. The penalty coefficient \code{lambda} is treated separately.
#' 
#' The log-likelihood functions \code{\link[evmix:fpsden]{lpsden}} and \code{\link[evmix:fpsden]{nlpsden}}
#' evaluate the likelihood for the original dataset, using the fitted P-splines density estimator. The
#' log-likelihood is output as \code{nllh} from the fitting function \code{\link[evmix:fpsden]{fpsden}}.
#' They do not provide the likelihood for the Poisson regression of the histogram counts, which is usually
#' evaluated using the deviance. The deviance (via CVMSE for Poisson counts) is also output as \code{cvlambda}
#' from the fitting function \code{\link[evmix:fpsden]{fpsden}}.
#' 
#' The \code{\link[evmix:fpsden]{iwlspsden}} function performs the IWLS. The 
#' \code{\link[evmix:fpsden]{cvpsden}} function calculates the leave-one-out cross-validation 
#' sum of the squared errors. They are not designed to be used directly by users. No checks of the
#' inputs are carried out.
#' 
#' @return Log-likelihood for original data is given by \code{\link[evmix:fpsden]{lpsden}} and it's
#'   wrappers for negative log-likelihood from \code{\link[evmix:fpsden]{nlpsden}}. Cross-validation 
#'   sum of square of errors is provided by \code{\link[evmix:fpsden]{cvpsden}}. Poisson regression
#'   fitting by IWLS is carried out in \code{\link[evmix:fpsden]{iwlspsden}}. Fitting function
#'   \code{\link[evmix:fpsden]{fpsden}} returns a simple list with the
#'   following elements
#'
#' \tabular{ll}{
#'  \code{call}:                \tab \code{optim} call\cr
#'  \code{x}:                   \tab data vector \code{x}\cr
#'  \code{xrange}:              \tab range of support of B-splines\cr
#'  \code{degree}:              \tab degree of B-splines\cr
#'  \code{nseg}:                \tab number of internal segments\cr
#'  \code{design.knots}:        \tab knots used in \code{\link[splines:splineDesign]{splineDesign}}\cr
#'  \code{ord}:                 \tab order of penalty term\cr
#'  \code{binned}:              \tab histogram results\cr
#'  \code{breaks}:              \tab histogram breaks\cr
#'  \code{mids}:                \tab histogram mid-bins\cr
#'  \code{counts}:              \tab histogram counts\cr
#'  \code{nbinwidth}:           \tab scaling factor to convert counts to density\cr
#'  \code{bsplines}:            \tab B-splines matrix used for binned counts\cr
#'  \code{databsplines}:        \tab B-splines matrix used for data\cr
#'  \code{counts}:              \tab histogram counts\cr
#'  \code{lambdaseq}:           \tab \eqn{\lambda} vector for profile likelihood or scalar for fixed \eqn{\lambda}\cr
#'  \code{cvlambda}:            \tab CV MSE for each \eqn{\lambda}\cr
#'  \code{mle} and \code{beta}: \tab vector of MLE of coefficients\cr
#'  \code{nllh}:                \tab negative log-likelihood for original data\cr
#'  \code{n}:                   \tab total original sample size\cr
#'  \code{lambda}:              \tab Estimated or fixed \eqn{\lambda}\cr
#' }
#' 
#' @note The data are both vectors. Infinite and missing sample values are dropped.
#' 
#' No initial values for the coefficients are needed.
#' 
#' It is advised to specify the range of support \code{xrange}, using finite end-points. This is 
#' especially important when the support is bounded. By default \code{xrange} is simply the range of the
#' input data \code{range(x)}.
#' 
#' Further, it is advised to always set the histogram bin \code{breaks}, expecially if the support is bounded.
#' By default \code{10*ln(n)} equi-spaced bins are defined between \code{xrange}.
#' 
#' @references
#' \url{http://www.math.canterbury.ac.nz/~c.scarrott/evmix}
#' 
#' \url{http://en.wikipedia.org/wiki/Cross-validation_(statistics)}
#' 
#' \url{http://en.wikipedia.org/wiki/B-spline}
#' 
#' \url{http://www.stat.lsu.edu/faculty/marx}
#' 
#' Eilers, P.H.C. and Marx, B.D. (1996). Flexible smoothing with B-splines and penalties.
#' Statistical Science 11(2), 89-121.
#' 
#' @author Alfadino Akbar and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @section Acknowledgments: The Poisson regression and leave-one-out cross-validation functions
#' are based on the code of Eilers and Marx (1996) available from Brian Marx's website 
#' \url{http://www.stat.lsu.edu/faculty/marx}, which is gratefully acknowledged.
#' 
#' @seealso \code{\link[evmix:kden]{kden}}.
#'  
#' @aliases fpsden lpsden nlpsden iwlspsden cvpsden
#' @family  psden fpsden
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' par(mfrow = c(1, 1))
#' 
#' x = rnorm(1000)
#' xx = seq(-4, 4, 0.01)
#' y = dnorm(xx)
#' 
#' # Plenty of histogram bins (100)
#' breaks = seq(-4, 4, length.out=101)
#' 
#' # P-spline fitting with cubic B-splines, 2nd order penalty and 10 internal segments
#' # CV search for penalty coefficient. 
#' fit = fpsden(x, lambdaseq = 10^seq(-5, 5, 0.25), breaks = breaks,
#'              xrange = c(-4, 4), nseg = 10, degree = 3, ord = 2)
#' psdensity = exp(fit$bsplines %*% fit$mle)
#' 
#' hist(x, freq = FALSE, breaks = seq(-4, 4, length.out=101), xlim = c(-6, 6))
#' lines(xx, y, col = "black") # true density
#' 
#' lines(fit$mids, psdensity/fit$nbinwidth, lwd = 2, col = "blue") # P-splines density
#' 
#' # check density against dpsden function
#' with(fit, lines(xx, dpsden(xx, beta, nbinwidth, design = design.knots),
#'                 lwd = 2, col = "red", lty = 2))
#'
#' # vertical lines for all knots
#' with(fit, abline(v = design.knots, col = "red"))
#'
#' # internal knots
#' with(fit, abline(v = design.knots[(degree + 2):(length(design.knots) - degree - 1)], col = "blue"))
#'   
#' # boundary knots (support of B-splines)
#' with(fit, abline(v = design.knots[c(degree + 1, length(design.knots) - degree)], col = "green"))
#'
#' legend("topright", c("True Density","P-spline density","Using dpsdens function"),
#'   col=c("black", "blue", "red"), lty = c(1, 1, 2))
#' legend("topleft", c("Internal Knots", "Boundaries", "Extra Knots"),
#'   col=c("blue", "green", "red"), lty = 1)
#' }
#'   

# maximum likelihood fitting for P-splines density estimate
fpsden <- function(x, lambdaseq = NULL, breaks = NULL, xrange = NULL,
                   nseg = 10, degree = 3, design.knots = NULL, ord = 2) {

  call <- match.call()
    
  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(lambdaseq, allowvec = TRUE, allownull = TRUE)
  
  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    x = x[is.finite(x)] # ignore missing and infinite cases
  }

  check.quant(x)
  n = length(x)

  # Check B-spline specification
  check.param(xrange, allowvec = TRUE, allownull = TRUE)
  check.n(nseg) # optional, if knots provided
  check.n(degree, allowzero = TRUE)
  check.param(design.knots, allowvec = TRUE, allownull = TRUE) # optional
  check.n(ord)

  # Two options for specifying knots:
  #    1) design.knots vector
  #    2) xrange, nseg and degree
  # if both provided then design.knots is used
  if (is.null(design.knots) & is.null(xrange)) {
    xrange = range(x)
    
    # wee buffer on support for machine precision issues in splineDesign range checks
    # But also checks for most common lower and upper bounds
    if (isTRUE(all.equal(xrange[1], 0)) & (xrange[1] > 0)) {
      xrange[1] = 0
    } else if (isTRUE(all.equal(xrange[1], -1)) & (xrange[1] > -1)) {
      xrange[1] = -1
    } else {
      xrange[1] = xrange[1] - 8*.Machine$double.eps    
    }
    if (isTRUE(all.equal(xrange[2], 1)) & (xrange[2] < 1)) {
      xrange[2] = 1
    } else if (isTRUE(all.equal(xrange[2], 100)) & (xrange[2] < 100)) {
      xrange[2] = 100
    } else {
      xrange[2] = xrange[2] + 8*.Machine$double.eps
    }
  }

  if (is.null(design.knots)) {
    # check x-range
    if (length(xrange) != 2) stop("knot range in xrange must be vector of length 2")
    if (diff(xrange) <= 0) stop("knot range in xrange must have positive width")
    
    # consistent with Eilers and Marx the "P-spline masters":
    # defaults to regular knots and not natural B-splines,
    # so each B-spline is just shifted/spliced version of each other
    dx = diff(xrange)/nseg # regular knot spacing
    design.knots = seq(xrange[1] - degree * dx, xrange[2] + degree * dx, by = dx)
    
  } else {
    # if knots specified, they must be sorted
    if (is.unsorted(design.knots)) {
      design.knots = sort(design.knots)
    } else {
      if (design.knots[1] > design.knots[length(design.knots)])
        design.knots = rev(design.knots)
    }
    
    # cannot check degree as beta coefficients not provided (as checked in psden functions)
    
    # number of segments and degree determine design.knots, and vice-versa, so stop if different
    # different behaviour to psden function, which can check degree in advance using beta
    if (nseg != (length(design.knots) - 1 - degree*2)) {
      stop(paste("Number of segments = ", nseg, " and degree = ", degree, 
                    " are inconsistent with design knots where length(design.knots) = ",
                    length(design.knots), " should be nseg + 1 + degree*2 = ",
                    ", so ignoring nseg and degree", sep=""))
    }

    # xrange also determined by design knots and degree
    if (!is.null(xrange)) {
      if (!isTRUE(all.equal(xrange, design.knots[c(degree + 1, length(design.knots) - degree)]))) {
        warning(paste("Interpolation range (", xrange[1], ", ", xrange[2], ") is inconsistent with design knots (",
                      design.knots[degree + 1], ", " , design.knots[length(design.knots) - degree],
                      "), so xrange reset", sep=""))
        xrange = design.knots[c(degree + 1, length(design.knots) - degree)]
      }
    } else {
      xrange = design.knots[c(degree + 1, length(design.knots) - degree)]
    }
  }
  
  # Calculate histogram counts
  if (is.null(breaks)) {
    breaks = ceiling(10 * log(n)) # Eilers and Marx default setting for # breaks
  } else {
    # if breaks specified, they must be sorted
    if (is.unsorted(breaks)) {
      breaks = sort(breaks)
    } else {
      if (breaks[1] > breaks[length(breaks)])
        breaks = rev(breaks)
    }
  }
  
  # If breaks are not prescribed, then set sequence from min to max of data
  # Different to Eilers and Marx (1996) who go slightly beyond the data
  # Our approach makes more sense when data is bounded (density cannot go beyond min and max of data)
  # Can be overridden by specifying xrange
  if (length(breaks) == 1) {
    breaks = seq(xrange[1], xrange[2], length.out = breaks)
  }

  # Slightly modified from Eilers and Marx (1996) code:
  #   1) x-locations are mid-points of bins rather than lower boundary
  #   2) if # of breaks is provided (rather than actual breaks) then 
  #      actual breaks are specified by default behaviour of hist function
    
  # breaks must cover all data (hist function will check this as well)
  if (any((x < breaks[1]) | (x > breaks[length(breaks)])))
    stop(paste("breaks range (", breaks[1], ", ", breaks[length(breaks)], 
               ") must cover range of input data", sep=""))

  binned = hist(x, breaks, plot = FALSE)
  mids = binned$mids
  counts = binned$counts
  binwidth = binned$breaks[2] - binned$breaks[1]
  nbinwidth = n * binwidth
  
  # B-spline basis at the bin mid-points
  bsplines = as.matrix.csr(splineDesign(design.knots, mids, ord = degree + 1))
  
  # Check if profile likelihood or fixed lambda is being used
  # and determine initial values for parameters in each case
  if (is.null(lambdaseq)) { # not profile or fixed
      lambda = 10 # same as Eilers and Marx (1996)
      cvlambda = cvpsden(lambda, counts = counts, bsplines = bsplines, ord = ord)
  } else { # profile or fixed
    
    # CVMSE for lambda or scalar given
    if (length(lambdaseq) != 1) {
      
      cvlambda = sapply(lambdaseq, cvpsden, counts = counts, bsplines = bsplines, ord = ord)
      
      if (all(!is.finite(cvlambda))) stop("lambdas are all invalid")
      lambda = lambdaseq[which.min(cvlambda)]

    } else {
      lambda = lambdaseq
      cvlambda = cvpsden(lambda, counts = counts, bsplines = bsplines, ord = ord)
    }
  }

  # IWLS for estimating B-spline coefficients, given penalty coefficient
  beta = iwlspsden(counts = counts, bsplines = bsplines, ord = ord, lambda = lambda)

  # B-spline basis at the actual data needed to evaluate log-likelihood
  databsplines = as.matrix.csr(splineDesign(design.knots, x, ord = degree + 1))

  nllh = nlpsden(beta, x, databsplines, nbinwidth)

  list(call = call, x = as.vector(x), xrange = xrange, degree = degree, nseg = nseg, 
       design.knots = design.knots, ord = ord, binned = binned,
       breaks = breaks, mids = mids, counts = counts, nbinwidth = nbinwidth, bsplines = bsplines,
       databsplines = databsplines, lambdaseq = lambdaseq, cvlambda = cvlambda,
       mle = beta, nllh = nllh, n = n, beta = beta, lambda = lambda)
}

#' @export
#' @aliases fpsden lpsden nlpsden iwlspsden cvpsden
#' @rdname  fpsden

# log-likelihood function for P-splines density estimate
lpsden <- function(x, beta = NULL, bsplines = NULL, nbinwidth = 1, log = TRUE) {

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(beta, allowvec = TRUE)
  check.posparam(nbinwidth, allowvec = TRUE)
  check.logic(log)
  
  if (!(is.matrix.csr(bsplines) | (is.matrix(bsplines))))
    stop("bsplines input must be matrix, with each column a B-spline")

  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    x = x[is.finite(x)] # ignore missing and infinite cases
  }

  check.quant(x)
  n = length(x)
  np = length(beta)
  
  if (dim(bsplines)[1] != n)
    stop("bsplines input must be matrix, with row for each datapoint in x")

  if (dim(bsplines)[2] != np)
    stop("bsplines input must be matrix, with column for each B-spline")

  pred = bsplines %*% beta
  l = sum(pred) - n * log(nbinwidth)
  
  if (!log) l = exp(l)  
  
  l
}

#' @export
#' @aliases fpsden lpsden nlpsden iwlspsden cvpsden
#' @rdname  fpsden
  
# negative log-likelihood function for P-splines density estimate
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlpsden <- function(pvector, x, bsplines = NULL, nbinwidth = 1, finitelik = FALSE) {

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(pvector, allowvec = TRUE)
  check.posparam(nbinwidth)
  check.logic(finitelik)

  if (!(is.matrix.csr(bsplines) | (is.matrix(bsplines))))
    stop("bsplines input must be matrix, with each column a B-spline")

  n = length(x)
  np = length(pvector)
    
  if (dim(bsplines)[1] != n)
    stop("bsplines input must be matrix, with row for each datapoint in x")

  if (dim(bsplines)[2] != np)
    stop("bsplines input must be matrix, with column for each B-spline")
  
  beta = pvector
  
  nllh = -lpsden(x, beta, bsplines, nbinwidth) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}

#' @export
#' @aliases fpsden lpsden nlpsden iwlspsden cvpsden
#' @rdname  fpsden

# leave one cross-validation RMSE for P-splines density estimate, with given B-splines
cvpsden = function(lambda = 1, counts, bsplines, ord = 2) {

  # Adapted from poisson.cv function of Paul Eilers and Brian Marx (2007)
  # Not designed to be called by users, so no input checking (to make it more efficient)
  
  nc = length(counts)
  np = dim(bsplines)[2]
  
  cv.pred <- function(i, counts, bsplines, ord, lambda) {
    exp(bsplines[i, ] %*% iwlspsden(counts[-i], bsplines[-i, ], ord, lambda))
  }

  if (lambda <= 0) {
    cv = Inf
  } else {
    pred = sapply(1:nc, cv.pred, counts = counts, bsplines = bsplines, ord = ord, lambda = lambda) 
    cv = mean((counts - pred)^2) # CVMSE
  }
  
  return(cv)
}

#' @export
#' @aliases fpsden lpsden nlpsden iwlspsden cvpsden
#' @rdname  fpsden

# iterative weighted least squares fitting for P-splines density estimate, with given B-splines
iwlspsden = function(counts, bsplines, ord = 2, lambda = 10) {

  # Adapted from ps.poisson function of Paul Eilers and Brian Marx (2007)
  # Not designed to be called by users, so no input checking (to make it more efficient)

  nc = length(counts)
  np = dim(bsplines)[2]

  # Penalty matrix from differences
  P = as.matrix.csr(sqrt(lambda) * diff(diag(np), diff = ord))
  X = rbind.matrix.csr(bsplines, P)

  npen = np - ord # same as dim(P)[1] (not checked)
  
  # Initialize
  newpred = log(counts + 0.0001)
  pred = rep(1, nc)

  # Iterations
  niter = 1
  diffpred = sum(abs(pred - newpred))
  while ((niter < 20) & (diffpred > 1e-5)) {
    mu = exp(pred)
    
    # determine weights and score
    w = ifelse(mu < 1e-16, 0, as.vector(mu)) # practical definition of zero predicted density

    newy = ifelse(w == 0, 0, as.vector((counts - mu) / w + pred))
    
    Y = c(newy, rep(0, npen))

    # Mixed model representation
    beta = slm.wfit(X, Y, weights = c(w, rep(1, npen)))$coef

    newpred = bsplines %*% beta # on log-link scale

    niter = niter + 1
    diffpred = sum(abs(pred - newpred))

    pred = newpred
  }
  
  if (diffpred > 1e-5) warning("IRLS fitting did not converge in 20 iterations, check all inputs")
  
  return(beta)
}
