#' @export
#' 
#' @title MLE Fitting of P-splines Density Estimate for Bulk and GPD Tail Extreme Value Mixture Model
#'
#' @description Maximum likelihood estimation for fitting the extreme value 
#' mixture model with P-splines density estimate for bulk distribution upto the threshold and conditional
#' GPD above threshold. With options for profile likelihood estimation for threshold and
#' fixed threshold approach.
#'
#' @param psdenx       P-splines based density estimate for each datapoint in x
#' @param bsplinefit   list output from P-splines density fitting \code{\link[evmix:fpsden]{fpsden}} function
#' @param phib         renormalisation constant for bulk model density \eqn{(1-\phi_u)/H(u)}, to make it integrate to \code{1-phiu}
#' @inheritParams fnormgpd
#' @inheritParams fgpd
#' @inheritParams fpsden
#' 
#' @details The extreme value mixture model with P-splines density estimate for bulk and GPD tail is 
#' fitted to the entire dataset. A two-stage maximum likelihood inference approach is taken. The first
#' stage consists fitting of the P-spline density estimator, which is acheived by MLE using the 
#' \code{\link[evmix:fpsden]{fpsden}} function. The second stage, conditions on the B-spline coefficients,
#' using MLE for the extreme value mixture model (GPD parameters and threshold, if requested). The estimated
#' parameters, variance-covariance matrix and their standard errors are automatically
#' output.
#' 
#' See help for \code{\link[evmix:fnormgpd]{fnormgpd}} for details of extreme value mixture models,
#' type \code{help fnormgpd}. Only the different features are outlined below for brevity.
#' 
#' As the second stage conditions on the Bs-pline coefficients, the full parameter vector is
#' (\code{u}, \code{sigmau}, \code{xi}) if threshold is also estimated and
#' (\code{sigmau}, \code{xi}) for profile likelihood or fixed threshold approach.
#' 
#' (Penalized) MLE estimation of the B-Spline coefficients is carried out using Poisson regression
#' based on histogram bin counts. See help for \code{\link[evmix:fpsden]{fpsden}} for details,
#' type \code{help fpsden}.
#' 
#' @return Log-likelihood is given by \code{\link[evmix:fpsdengpd]{lpsdengpd}} and it's
#'   wrappers for negative log-likelihood from \code{\link[evmix:fpsdengpd]{nlpsdengpd}}
#'   and \code{\link[evmix:fpsdengpd]{nlupsdengpd}}. Profile likelihood for single
#'   threshold given by \code{\link[evmix:fpsdengpd]{proflupsdengpd}}. Fitting function
#'   \code{\link[evmix:fpsdengpd]{fpsdengpd}} returns a simple list with the
#'   following elements
#'
#' \tabular{ll}{
#'  \code{call}:          \tab \code{optim} call\cr
#'  \code{x}:             \tab data vector \code{x}\cr
#'  \code{init}:          \tab \code{pvector}\cr
#'  \code{fixedu}:        \tab fixed threshold, logical\cr
#'  \code{useq}:          \tab threshold vector for profile likelihood or scalar for fixed threshold\cr
#'  \code{nllhuseq}:      \tab profile negative log-likelihood at each threshold in useq\cr
#'  \code{bsplinefit}:    \tab complete \code{\link[evmix:fpsden]{fpsden}} output\cr
#'  \code{psdenx}:        \tab P-splines based density estimate for each datapoint in \code{x}\cr
#'  \code{xrange}:        \tab range of support of B-splines\cr
#'  \code{degree}:        \tab degree of B-splines\cr
#'  \code{nseg}:          \tab number of internal segments\cr
#'  \code{design.knots}:  \tab knots used in \code{\link[splines:splineDesign]{splineDesign}}\cr
#'  \code{nbinwidth}:     \tab scaling factor to convert counts to density\cr
#'  \code{optim}:         \tab complete \code{optim} output\cr
#'  \code{conv}:          \tab indicator for "possible" convergence\cr
#'  \code{mle}:           \tab vector of MLE of (GPD and threshold, if relevant) parameters\cr
#'  \code{cov}:           \tab variance-covariance matrix of MLE of parameters\cr
#'  \code{se}:            \tab vector of standard errors of MLE of parameters\cr
#'  \code{rate}:          \tab \code{phiu} to be consistent with \code{\link[evd:fpot]{evd}}\cr
#'  \code{nllh}:          \tab minimum negative log-likelihood\cr
#'  \code{n}:             \tab total sample size\cr
#'  \code{beta}:          \tab vector of MLE of B-spline coefficients\cr
#'  \code{lambda}:        \tab Estimated or fixed \eqn{\lambda}\cr
#'  \code{u}:             \tab threshold (fixed or MLE)\cr
#'  \code{sigmau}:        \tab MLE of GPD scale\cr
#'  \code{xi}:            \tab MLE of GPD shape\cr
#'  \code{phiu}:          \tab MLE of tail fraction (bulk model or parameterised approach)\cr
#'  \code{se.phiu}:       \tab standard error of MLE of tail fraction\cr
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
#' When \code{pvector=NULL} then the initial values are:
#' \itemize{
#'  \item threshold 90\% quantile (not relevant for profile likelihood for threshold or fixed threshold approaches);
#'  \item MLE of GPD parameters above threshold. 
#' }
#' 
#' @references
#' \url{http://www.math.canterbury.ac.nz/~c.scarrott/evmix}
#' 
#' \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
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
#' Scarrott, C.J. and MacDonald, A. (2012). A review of extreme value
#' threshold estimation and uncertainty quantification. REVSTAT - Statistical
#' Journal 10(1), 33-59. Available from \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}
#' 
#' @author Alfadino Akbar and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @section Acknowledgments: See Acknowledgments in
#'   \code{\link[evmix:fnormgpd]{fnormgpd}}, type \code{help fnormgpd}.
#'   
#'   The Poisson regression and leave-one-out cross-validation functions
#' are based on the code of Eilers and Marx (1996) available from Brian Marx's website 
#' \url{http://www.stat.lsu.edu/faculty/marx}, which is gratefully acknowledged.
#' 
#' @seealso \code{\link[evmix:fpsden]{fpsden}}, \code{\link[evmix:normgpd]{fnormgpd}},
#' \code{\link[evmix:fgpd]{fgpd}} and \code{\link[evmix:gpd]{gpd}}
#'  
#' @aliases fpsdengpd lpsdengpd nlpsdengpd proflupsdengpd nlupsdengpd
#' @family  psdengpd fpsdengpd normgpd fnormgpd psden fpsden
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
#' fit = fpsdengpd(x, useq = seq(0, 3, 0.1), fixedu = TRUE,
#'              lambdaseq = 10^seq(-5, 5, 0.25), breaks = breaks,
#'              xrange = c(-4, 4), nseg = 10, degree = 3, ord = 2)
#'              
#' hist(x, freq = FALSE, breaks = breaks, xlim = c(-6, 6))
#' lines(xx, y, col = "black") # true density
#' 
#' # P-splines+GPD
#' with(fit, lines(xx, dpsdengpd(xx, beta, nbinwidth, 
#'                               u = u, sigmau = sigmau, xi = xi, design = design.knots),
#'                 lwd = 2, col = "red"))
#' abline(v = fit$u, col = "red", lwd = 2, lty = 3)
#' 
#' # P-splines density estimate
#' with(fit, lines(xx, dpsden(xx, beta, nbinwidth, design = design.knots),
#'                 lwd = 2, col = "blue", lty = 2))
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
#' legend("topright", c("True Density","P-spline density","P-spline+GPD"),
#'   col=c("black", "blue", "red"), lty = c(1, 2, 1))
#' legend("topleft", c("Internal Knots", "Boundaries", "Extra Knots", "Threshold"),
#'   col=c("blue", "green", "red", "red"), lty = c(1, 1, 1, 2))
#' }
#'   

# maximum likelihood fitting for P-splines density estimate for bulk with GPD for upper tail
fpsdengpd <- function(x, phiu = TRUE, useq = NULL, fixedu = FALSE, pvector = NULL, 
                      lambdaseq = NULL, breaks = NULL, xrange = NULL,
                      nseg = 10, degree = 3, design.knots = NULL, ord = 2, 
                      std.err = TRUE, method = "BFGS", control = list(maxit = 10000), finitelik = TRUE, ...) {

  call <- match.call()
    
  np = 3 # maximum number of parameters

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.logic(phiu)
  check.param(useq, allowvec = TRUE, allownull = TRUE)
  check.logic(fixedu)
  check.logic(std.err)
  check.optim(method)
  check.control(control)
  check.logic(finitelik)

  check.param(lambdaseq, allowvec = TRUE, allownull = TRUE)
  
  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    x = x[is.finite(x)] # ignore missing and infinite cases
  }

  check.quant(x)
  n = length(x)

  if ((method == "L-BFGS-B") | (method == "BFGS")) finitelik = TRUE

  # Likelihood estimation is done "conditional" on P-spline fit (two stage approach)
  
  # Stage 1 - Use fpsden to do P-spline fitting first
  
  bsplinefit = fpsden(x, lambdaseq = lambdaseq, breaks = breaks, xrange = xrange,
         nseg = nseg, degree = degree, design.knots = design.knots, ord = ord)
  
  # P-spline based density estimate at each data point
  psdenx = with(bsplinefit, exp(databsplines %*% beta) / nbinwidth) 

  # Stage 2 - extreme value modelling conditional on B-splines coefficients
  
  # useq must be specified if threshold is fixed
  if (fixedu & is.null(useq))
    stop("for fixed threshold approach, useq must be specified (as scalar or vector)")

  # Check if profile likelihood or fixed threshold is being used
  # and determine initial values for parameters in each case
  if (is.null(useq)) { # not profile or fixed

    check.nparam(pvector, nparam = np, allownull = TRUE)
    
    if (is.null(pvector)) {
      pvector[1] = as.vector(quantile(x, 0.9))
      initfgpd = fgpd(x, pvector[1], std.err = FALSE)
      pvector[2] = initfgpd$sigmau
      pvector[3] = initfgpd$xi
    }
    
  } else { # profile or fixed
    
    check.nparam(pvector, nparam = np - 1, allownull = TRUE)

    # profile likelihood for threshold or scalar given
    if (length(useq) != 1) {
      
      # remove thresholds with less than 5 excesses
      useq = useq[sapply(useq, FUN = function(u, x) sum(x > u) > 5, x = x)]
      check.param(useq, allowvec = TRUE)
      
      nllhu = sapply(useq, proflupsdengpd, pvector = pvector, x = x, 
                     psdenx = psdenx, phiu = phiu, bsplinefit = bsplinefit,
                     method = method, control = control, finitelik = finitelik, ...) # do not need phib
      
      if (all(!is.finite(nllhu))) stop("thresholds are all invalid")
      u = useq[which.min(nllhu)]

      if ((which.min(nllhu) == 1) | (which.min(nllhu) == length(nllhu))) 
        warning("Profile MLE of threshold is on boundary of sequence (useq)")
    } else {
      u = useq
    }

    if (fixedu) { # threshold fixed
      if (is.null(pvector)) {
        initfgpd = fgpd(x, u, std.err = FALSE)
        pvector[1] = initfgpd$sigmau
        pvector[2] = initfgpd$xi
      }
    } else { # threshold as initial value in usual MLE
      if (is.null(pvector)) {
        pvector[1] = u
        initfgpd = fgpd(x, pvector[1], std.err = FALSE)
        pvector[2] = initfgpd$sigmau
        pvector[3] = initfgpd$xi
      } else {
        pvector = c(u, pvector) # prefix with estimated threshold
      }
    }
  }

  if (fixedu) { # fixed threshold (separable) likelihood

    # For given u, phiu can be determined as we are conditioning on P-splines DE
    #     phiu == TRUE it comes from bulk model
    #     phiu == FALSE it comes from new parameter
    # Usually phiu is updated for bulk nmodel parameters, but due to conditioning
    # we can pre-define it in the likelihood
    # Gives computational efficiency as avoids integrating P-splines DE upto u on
    # each optimisation step
    if (is.logical(phiu)) {
      pu = with(bsplinefit, try(integrate(pscounts, lower = u, upper = xrange[2],
                         beta = beta, design.knots = design.knots, degree = degree,
                         subdivisions = 10000, rel.tol = 1e-9, stop.on.error = FALSE)))
    
      if (inherits(pu, "try-error")) {
        pu$value = NA
        stop("failed to numerically evaluate cdf of P-spline density estimate")
      }
      pu = 1 - pu$value/bsplinefit$nbinwidth
    
      if (phiu) {
        phiu = 1 - pu
        se.phiu = NA
      } else {
        phiu = mean(x > u, na.rm = TRUE) # infinite cases already dropped
        se.phiu = sqrt(phiu * (1 - phiu) / n)
      }
    }
    phib = (1 - phiu) / pu
    
    nllh = nlupsdengpd(pvector, u, x, psdenx, phiu, bsplinefit, phib)
    if (is.infinite(nllh)) {
      pvector[2] = 0.1
      nllh = nlupsdengpd(pvector, u, x, psdenx, phiu, bsplinefit, phib)
    }
    if (is.infinite(nllh)) stop("initial parameter values are invalid")
  
    fit = optim(par = as.vector(pvector), fn = nlupsdengpd,
                u = u, x = x, psdenx = psdenx, phiu = phiu, bsplinefit = bsplinefit, phib = phib,
                finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)

    sigmau = fit$par[1]
    xi = fit$par[2]
    
  } else { # complete (non-separable) likelihood
    
    nllh = nlpsdengpd(pvector, x, psdenx, phiu, bsplinefit = bsplinefit) # do not need phib
    if (is.infinite(nllh)) {
      pvector[3] = 0.1
      nllh = nlpsdengpd(pvector, x, psdenx, phiu, bsplinefit, phib)
    }
    if (is.infinite(nllh)) stop("initial parameter values are invalid")
  
    fit = optim(par = as.vector(pvector), fn = nlpsdengpd,
                x = x, psdenx = psdenx, phiu = phiu, bsplinefit = bsplinefit,
                finitelik = finitelik, method = method, control = control, hessian = TRUE, ...) # do not need phib
    
    u = fit$par[1]
    sigmau = fit$par[2]
    xi = fit$par[3]

    pu = with(bsplinefit, try(integrate(pscounts, lower = u, upper = xrange[2],
      beta = beta, design.knots = design.knots, degree = degree,
      subdivisions = 10000, rel.tol = 1e-9, stop.on.error = FALSE)))
    
    if (inherits(pu, "try-error")) {
      pu$value = NA
      stop("failed to numerically evaluate cdf of P-spline density estimate")
    }
    pu = 1 - pu$value/bsplinefit$nbinwidth
    
    if (phiu) {
      phiu = 1 - pu
      se.phiu = NA
    } else {
      phiu = mean(x > u, na.rm = TRUE)
      se.phiu = sqrt(phiu * (1 - phiu) / n)
    }
  }

  conv = TRUE
  if ((fit$convergence != 0) | any(fit$par == pvector) | (abs(fit$value) >= 1e6)) {
    conv = FALSE
    warning("check convergence")
  }

  if (std.err) {
    qrhess = qr(fit$hessian)
    if (qrhess$rank != ncol(qrhess$qr)) {
      warning("observed information matrix is singular")
      se = NULL
      invhess = NULL
    } else {
      invhess = solve(qrhess)
      vars = diag(invhess)
      if (any(vars <= 0)) {
        warning("observed information matrix is singular")
        invhess = NULL
        se = NULL
      } else {
        se = sqrt(vars)
      }  
    }
  } else {
    invhess = NULL
    se = NULL
  }
  
  if (!exists("nllhu")) nllhu = NULL
    
  list(call = call, x = as.vector(x), init = as.vector(pvector), fixedu = fixedu, useq = useq, nllhuseq = nllhu,
       bsplinefit = bsplinefit, psdenx = psdenx,
       xrange = bsplinefit$xrange, degree = bsplinefit$degree, nseg = bsplinefit$nseg, 
       design.knots = bsplinefit$design.knots, nbinwidth = bsplinefit$nbinwidth,      
       optim = fit, conv = conv, cov = invhess, mle = fit$par, se = se, rate = phiu,
       nllh = fit$value, n = n, beta = bsplinefit$beta, lambda = bsplinefit$lambda,
       u = u, sigmau = sigmau, xi = xi, phiu = phiu, se.phiu = se.phiu)
}

#' @export
#' @aliases fpsdengpd lpsdengpd nlpsdengpd proflupsdengpd nlupsdengpd
#' @rdname  fpsdengpd

# log-likelihood function for "conditional" P-splines density estimate for bulk with GPD for upper tail
lpsdengpd <- function(x, psdenx, u = NULL, sigmau = NULL, xi = 0, phiu = TRUE, bsplinefit = NULL,
                      phib = NULL, log = TRUE) {
  
  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(psdenx, allowvec = TRUE)
  check.param(u, allownull = TRUE) # defaults below
  check.param(sigmau, allownull = TRUE) # defaults below
  check.param(xi)
  check.phiu(phiu, allowfalse = TRUE)
  check.posparam(phib, allownull = TRUE) # fixed u case, see fitting function
  check.logic(log)
  
  if (missing(bsplinefit)) 
    stop("bsplinefit input must be given. It is the list that results from applying fpsden to your data.")
    
  if (length(psdenx) != length(x))
    stop("P-spline density estimate needed for all x, so length(x) == length(psdenx)")

  if ((!is.null(phib)) & (is.logical(phib)))
    stop("phib is numeric, so phiu must also be numeric (currently it is logical)")
    
  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    psdenx = psdenx[is.finite(x)] # ignore missing and infinite cases
    x = x[is.finite(x)] # ignore missing and infinite cases
  }

  # default for u and sigmau
  if (is.null(u)) {
    u = with(bsplinefit, qpsden(0.9, beta, nbinwidth, xrange, nseg, degree, design.knots))
    check.param(u)
  }
  if (is.null(sigmau)) {
    sigmau = diff(with(bsplinefit, qpsden(c(0.975, 0.85), beta, nbinwidth, xrange, nseg, degree,
                                          design.knots))) # if normal then approx one sd
    check.posparam(sigmau)
  }

  check.quant(x)
  n = length(x)
  np = length(bsplinefit$beta)
  
  check.inputn(c(length(u), length(sigmau), length(xi), length(phiu)), allowscalar = TRUE)
  if (!is.null(phib)) check.inputn(c(length(phib), length(phiu)), allowscalar = TRUE)

  # assume NA or NaN are irrelevant as entire lower tail is now modelled
  # inconsistent with evd library definition
  # hence use which() to ignore these

  whichu = which(x > u)
  xu = x[whichu]
  nu = length(xu)
  whichb = which(x <= u)
  xb = x[whichb]
  nb = length(xb)

  if (n != nb + nu) {
    stop("total non-finite sample size is not equal to those above threshold and those below or equal to it")    
  }

  if ((sigmau <= 0) | (u <= min(x)) | (u >= max(x))) {
    l = -Inf
  } else {
    if (is.null(phib)) {# usual (not fixed u) case
      if (is.logical(phiu)) {
        pu = with(bsplinefit, try(integrate(pscounts, lower = u, upper = xrange[2],
                                            beta = beta, design.knots = design.knots, degree = degree,
                                            subdivisions = 10000, rel.tol = 1e-9, stop.on.error = FALSE)))
    
        if (inherits(pu, "try-error")) {
          pu$value = NA
          stop("failed to numerically evaluate cdf of P-spline density estimate")
        }
        pu = 1 - pu$value/bsplinefit$nbinwidth
    
        if (phiu) {
          phiu = 1 - pu
        } else {
          phiu = nu / n
        }
      }
      phib = (1 - phiu) / pu
    } else { # fixed u case, where phib is predefined in fitting function
      pu = (1 - phiu) / phib
    }

    syu = 1 + xi * (xu - u) / sigmau  
  
    if ((min(syu) <= 0) | (phiu <= 0) | (phiu >= 1) | (pu <= 0) | (pu >= 1)) {
      l = -Inf
    } else { 
      l = lgpd(xu, u, sigmau, xi, phiu)
      l = l + sum(log(psdenx[whichb])) + nb*log(phib)
    }
  }
  
  if (!log) l = exp(l)
  
  l
}

#' @export
#' @aliases fpsdengpd lpsdengpd nlpsdengpd proflupsdengpd nlupsdengpd
#' @rdname  fpsdengpd

# negative log-likelihood function for P-splines density estimate for bulk with GPD for upper tail
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlpsdengpd <- function(pvector, x, psdenx, phiu = TRUE, bsplinefit, phib = NULL, finitelik = FALSE) {
  
  np = 3 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.phiu(phiu, allowfalse = TRUE)
  check.posparam(phib, allownull = TRUE)
  check.logic(finitelik)

  if (missing(bsplinefit)) 
    stop("bsplinefit input must be given. It is the list that results from applying fpsden to your data.")

  u = pvector[1]
  sigmau = pvector[2]
  xi = pvector[3]

  nllh = -lpsdengpd(x, psdenx, u, sigmau, xi, phiu, bsplinefit, phib)
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}

#' @export
#' @aliases fpsdengpd lpsdengpd nlpsdengpd proflupsdengpd nlupsdengpd
#' @rdname  fpsdengpd

# profile negative log-likelihood function for given threshold for
# P-splines density estimate for bulk with GPD for upper tail
# designed for sapply to loop over vector of thresholds (hence u is first input)
proflupsdengpd <- function(u, pvector, x, psdenx, phiu = TRUE, bsplinefit,
                          method = "BFGS", control = list(maxit = 10000), finitelik = TRUE, ...) {

  np = 3 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np - 1, allownull = TRUE)
  check.param(u)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(psdenx, allowvec = TRUE)
  check.phiu(phiu, allowfalse = TRUE)
  check.optim(method)
  check.control(control)
  check.logic(finitelik)

  if (missing(bsplinefit)) 
    stop("bsplinefit input must be given. It is the list that results from applying fpsden to your data.")  
  
  if (length(psdenx) != length(x))
    stop("P-spline density estimate needed for all x, so length(x) == length(psdenx)")

  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    x = x[is.finite(x)] # ignore missing and infinite cases
  }

  check.quant(x)
  
  # For given u, phiu can be determined as we are conditioning on P-splines DE
  # see fpsden for details
  if (is.logical(phiu)) {
    pu = with(bsplinefit, try(integrate(pscounts, lower = u, upper = xrange[2],
                                        beta = beta, design.knots = design.knots, degree = degree,
                                        subdivisions = 10000, rel.tol = 1e-9, stop.on.error = FALSE)))

    if (inherits(pu, "try-error")) {
      pu$value = NA
      stop("failed to numerically evaluate cdf of P-spline density estimate")
    }
    pu = 1 - pu$value/bsplinefit$nbinwidth
    
    if (phiu) {
      phiu = 1 - pu
    } else {
      phiu = mean(x > u, na.rm = TRUE) # infinite cases already dropped
    }
  }
  phib = (1 - phiu) / pu

  # check initial values for other parameters, try usual alternative
  if (!is.null(pvector)) {
    nllh = nlupsdengpd(pvector, u, x, psdenx, phiu, bsplinefit, phib)
    
    if (is.infinite(nllh)) pvector = NULL
  }

  if (is.null(pvector)) {
    initfgpd = fgpd(x, u, std.err = FALSE)
    pvector[1] = initfgpd$sigmau
    pvector[2] = initfgpd$xi
    nllh = nlupsdengpd(pvector, u, x, psdenx, phiu, bsplinefit, phib)
  }

  if (is.infinite(nllh)) {
    pvector[2] = 0.1
    nllh = nlupsdengpd(pvector, u, x, psdenx, phiu, bsplinefit, phib)
  }

  # if still invalid then output cleanly
  if (is.infinite(nllh)) {
    warning(paste("initial parameter values for threshold u =", u, "are invalid"))
    fit = list(par = rep(NA, np), value = Inf, counts = 0, convergence = NA, 
      message = "initial values invalid", hessian = rep(NA, np))
  } else {
    fit = optim(par = as.vector(pvector), fn = nlupsdengpd, u = u, x = x, psdenx = psdenx, phiu = phiu,
                       bsplinefit = bsplinefit, phib = phib,
                       finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)
  }
    
  if (finitelik & is.infinite(fit$value)) {
    fit$value = sign(fit$value) * 1e6
  }

  fit$value
}

#' @export
#' @aliases fpsdengpd lpsdengpd nlpsdengpd proflupsdengpd nlupsdengpd
#' @rdname  fpsdengpd

# negative log-likelihood function for P-splines density estimate for bulk with GPD for upper tail
# (wrapper for likelihood, designed for threshold to be fixed and other parameters optimised)
nlupsdengpd <- function(pvector, u, x, psdenx, phiu = TRUE, bsplinefit = bsplinefit, phib = NULL, finitelik = FALSE) {

  np = 3 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np - 1)
  check.param(u)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(psdenx, allowvec = TRUE)
  check.phiu(phiu, allowfalse = TRUE)
  check.posparam(phib, allownull = TRUE)
  check.logic(finitelik)

  if (missing(bsplinefit))
    stop("bsplinefit input must be given. It is the list that results from applying fpsden to your data.")
  
  sigmau = pvector[1]
  xi = pvector[2]

  nllh = -lpsdengpd(x, psdenx, u, sigmau, xi, phiu, bsplinefit, phib)
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}
