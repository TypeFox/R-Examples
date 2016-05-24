#' @export
#' 
#' @title MLE Fitting of Normal Bulk and GPD for Both Tails Interval Transition Mixture Model
#'
#' @description Maximum likelihood estimation for fitting the extreme value 
#' mixture model with normal for bulk distribution between thresholds, conditional
#' GPDs beyond thresholds and interval transition. With options for profile likelihood
#' estimation for both thresholds and interval half-width, which can also be fixed.
#'
#' @param ulseq    vector of lower thresholds (or scalar) to be considered in profile likelihood or
#'                \code{NULL} for no profile likelihood
#' @param urseq    vector of upper thresholds (or scalar) to be considered in profile likelihood or
#'                \code{NULL} for no profile likelihood
#' @param eulr    vector of epsilon, lower and upper thresholds considered in profile likelihood
#' @inheritParams fitmnormgpd
#' @inheritParams fitmgng
#' @inheritParams ditmgng
#' @inheritParams fnormgpd
#' @inheritParams fgpd
#' 
#' @details The extreme value mixture model with the normal bulk and GPD for both tails interval
#' transition is fitted to the entire dataset using maximum likelihood estimation. The estimated
#' parameters, variance-covariance matrix and their standard errors are automatically
#' output.
#' 
#' See \code{\link[evmix:itmgng]{ditmgng}} for explanation of GPD-normal-GPD interval
#' transition model, including mixing functions.
#' 
#' See also help for \code{\link[evmix:fnormgpd]{fnormgpd}} for details, type \code{help fnormgpd}. 
#' Only the different features are outlined below for brevity.
#' 
#' The full parameter vector is
#' (\code{nmean}, \code{nsd}, \code{epsilon}, \code{ul}, \code{sigmaul}, \code{xil},
#'                                            \code{ur}, \code{sigmaur}, \code{xir})
#' if thresholds and interval half-width are also estimated and
#' (\code{nmean}, \code{nsd}, \code{sigmaul}, \code{xil}, \code{sigmaur}, \code{xir})
#' for profile likelihood or fixed threshold approach.
#' 
#' If the profile likelihood approach is used, then a grid search over all combinations of epsilons and both thresholds
#' are carried out. The combinations which lead to less than 5 in any component outside of the
#' intervals are not considered.
#' 
#' A fixed pair of thresholds and epsilon approach is acheived by setting a single
#' scalar value to each in \code{ulseq}, \code{urseq} and \code{eseq} respectively.
#' 
#' @return Log-likelihood is given by \code{\link[evmix:fitmgng]{litmgng}} and it's
#'   wrappers for negative log-likelihood from \code{\link[evmix:fitmgng]{nlitmgng}}
#'   and \code{\link[evmix:fitmgng]{nluitmgng}}. Profile likelihood for 
#'   thresholds and interval half-width given by \code{\link[evmix:fitmgng]{profluitmgng}}.
#'   Fitting function \code{\link[evmix:fitmgng]{fitmgng}} returns a simple list with the
#'   following elements
#'
#' \tabular{ll}{
#'  \code{call}:      \tab \code{optim} call\cr
#'  \code{x}:         \tab data vector \code{x}\cr
#'  \code{init}:      \tab \code{pvector}\cr
#'  \code{fixedeu}:   \tab fixed epsilon and threshold, logical\cr
#'  \code{ulseq}:     \tab lower threshold vector for profile likelihood or scalar for fixed threshold\cr
#'  \code{urseq}:     \tab upper threshold vector for profile likelihood or scalar for fixed threshold\cr
#'  \code{eseq}:      \tab interval half-width vector for profile likelihood or scalar for fixed threshold\cr
#'  \code{nllheuseq}: \tab profile negative log-likelihood at each combination in (eseq, ulseq, urseq)\cr
#'  \code{optim}:     \tab complete \code{optim} output\cr
#'  \code{mle}:       \tab vector of MLE of parameters\cr
#'  \code{cov}:       \tab variance-covariance matrix of MLE of parameters\cr
#'  \code{se}:        \tab vector of standard errors of MLE of parameters\cr
#'  \code{nllh}:      \tab minimum negative log-likelihood\cr
#'  \code{n}:         \tab total sample size\cr
#'  \code{nmean}:     \tab MLE of normal mean\cr
#'  \code{nsd}:       \tab MLE of normal standard deviation\cr
#'  \code{epsilon}:   \tab MLE of transition half-width\cr
#'  \code{ul}:        \tab lower threshold (fixed or MLE)\cr
#'  \code{sigmaul}:   \tab MLE of lower tail GPD scale\cr
#'  \code{xil}:       \tab MLE of lower tail GPD shape\cr
#'  \code{ur}:        \tab upper threshold (fixed or MLE)\cr
#'  \code{sigmaur}:   \tab MLE of upper tail GPD scale\cr
#'  \code{xir}:       \tab MLE of upper tail GPD shape\cr
#' }
#' 
#' @note When \code{pvector=NULL} then the initial values are:
#' \itemize{
#'  \item MLE of normal parameters assuming entire population is normal; and
#'  \item lower threshold 10\% quantile (not relevant for profile likelihood for threshold or fixed threshold approaches);
#'  \item upper threshold 90\% quantile (not relevant for profile likelihood for threshold or fixed threshold approaches);
#'  \item MLE of GPD parameters beyond threshold. 
#' }
#' 
#' @references
#' \url{http://www.math.canterbury.ac.nz/~c.scarrott/evmix}
#' 
#' \url{http://en.wikipedia.org/wiki/Normal_distribution}
#' 
#' \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
#' 
#' Scarrott, C.J. and MacDonald, A. (2012). A review of extreme value
#' threshold estimation and uncertainty quantification. REVSTAT - Statistical
#' Journal 10(1), 33-59. Available from \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}
#' 
#' Holden, L. and Haug, O. (2013). A mixture model for unsupervised tail
#' estimation. arxiv:0902.4137
#' 
#' @author Alfadino Akbar and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @section Acknowledgments: See Acknowledgments in
#'   \code{\link[evmix:fnormgpd]{fnormgpd}}, type \code{help fnormgpd}. Based on code
#' by Xin Zhao produced for MATLAB.
#' 
#' @seealso \code{\link[evmix:fgng]{fgng}}, \code{\link[stats:Normal]{dnorm}},
#'  \code{\link[evmix:fgpd]{fgpd}} and \code{\link[evmix:gpd]{gpd}}
#'  
#' @aliases fitmgng litmgng nlitmgng profluitmgng nluitmgng
#' @family  normgpd normgpdcon gng gngcon fnormgpd fnormgpdcon fgng fgngcon
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
#' # MLE for complete parameter set (not recommended!)
#' fit = fitmgng(x)
#' hist(x, breaks = seq(-6, 6, 0.1), freq = FALSE, xlim = c(-4, 4))
#' lines(xx, y)
#' with(fit, lines(xx, ditmgng(xx, nmean, nsd, epsilon, ul, sigmaul, xil,
#'                                                      ur, sigmaur, xir), col="red"))
#' abline(v = fit$ul + fit$epsilon * seq(-1, 1), col = "red")
#' abline(v = fit$ur + fit$epsilon * seq(-1, 1), col = "darkred")
#'   
#' # Profile likelihood for threshold which is then fixed
#' fitfix = fitmgng(x, eseq = seq(0, 2, 0.1), ulseq = seq(-2.5, 0, 0.25), 
#'                                          urseq = seq(0, 2.5, 0.25), fixedeu = TRUE)
#' with(fitfix, lines(xx, ditmgng(xx, nmean, nsd, epsilon, ul, sigmaul, xil,
#'                                                       ur, sigmaur, xir), col="blue"))
#' abline(v = fitfix$ul + fitfix$epsilon * seq(-1, 1), col = "blue")
#' abline(v = fitfix$ur + fitfix$epsilon * seq(-1, 1), col = "darkblue")
#' legend("topright", c("True Density", "GPD-normal-GPD ITM", "Profile likelihood"),
#'   col=c("black", "red", "blue"), lty = 1)
#' }
#'   

# maximum likelihood fitting for normal bulk with GPD for both tails
# interval transition mixture model
fitmgng <- function(x, eseq = NULL, ulseq = NULL, urseq = NULL, fixedeu = FALSE, pvector = NULL,
  std.err = TRUE, method = "BFGS", control = list(maxit = 10000), finitelik = TRUE, ...) {

  call <- match.call()
    
  np = 9 # maximum number of parameters

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.posparam(eseq, allowvec = TRUE, allownull = TRUE, allowzero = TRUE)
  check.param(ulseq, allowvec = TRUE, allownull = TRUE)
  check.param(urseq, allowvec = TRUE, allownull = TRUE)
  check.logic(fixedeu)
  check.logic(std.err)
  check.optim(method)
  check.control(control)
  check.logic(finitelik)

  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    x = x[is.finite(x)] # ignore missing and infinite cases
  }

  check.quant(x)
  n = length(x)

  if ((method == "L-BFGS-B") | (method == "BFGS")) finitelik = TRUE
  
  # useq must be specified if epsilon is fixed
  if (fixedeu & is.null(eseq))
    stop("for fixed epsilon approach, eseq must be specified (as scalar or vector)")

  # useq must be specified if threshold is fixed
  if (fixedeu & (is.null(ulseq) | is.null(urseq)))
    stop("for fixed threshold approach, ulseq and urseq must be specified (as scalar or vector)")

  # Check if profile likelihood or fixed threshold is being used
  # and determine initial values for parameters in each case
  if (is.null(eseq) | is.null(ulseq) | is.null(urseq)) { # not profile or fixed

    check.nparam(pvector, nparam = np, allownull = TRUE)
    
    if (is.null(pvector)) {
      pvector[1] = mean(x, trim = 0.2)
      pvector[2] = sd(x)
      pvector[3] = pvector[2]
      pvector[4] = as.vector(quantile(x, 0.1))
      initfgpd = fgpd(-x, -pvector[4], std.err = FALSE)
      pvector[5] = initfgpd$sigmau
      pvector[6] = initfgpd$xi
      pvector[7] = as.vector(quantile(x, 0.9))
      initfgpd = fgpd(x, pvector[7], std.err = FALSE)
      pvector[8] = initfgpd$sigmau
      pvector[9] = initfgpd$xi
    }

  } else { # profile or fixed
    
    check.nparam(pvector, nparam = np - 3, allownull = TRUE)

    # profile likelihood for threshold or scalar given
    if ((length(eseq) != 1) | (length(ulseq) != 1) | (length(urseq) != 1)) {
      
      eugrid = expand.grid(eseq, ulseq, urseq)
      
      # remove combinations where interval is beyond range of data (must be at least 5 below lower interval)
      eugrid = eugrid[sapply(eugrid[, 2] - eugrid[, 1], FUN = function(u, x) sum(x < u) >= 5, x = x),]
      
      # remove combinations where interval is beyond range of data (must be at least 5 above upper interval)
      eugrid = eugrid[sapply(eugrid[, 3] + eugrid[, 1], FUN = function(u, x) sum(x > u) >= 5, x = x),]

      # remove combinations where bulk has too little data (must be at least 5 in between intervals)
      eugrid = eugrid[apply(cbind(eugrid[, 2] + eugrid[, 1], eugrid[, 3] - eugrid[, 1]), 1,
                             FUN = function(interval, x) sum((x > interval[1]) & (x < interval[2])) >= 5, x = x),]

      check.posparam(eugrid[, 1], allowvec = TRUE, allowzero = TRUE)
      check.param(eugrid[, 2], allowvec = TRUE)
      check.param(eugrid[, 3], allowvec = TRUE)
      
      nllheu = apply(eugrid, 1, profleuitmgng, pvector = pvector, x = x,
                     method = method, control = control, finitelik = finitelik, ...)
      
      if (all(!is.finite(nllheu))) stop("thresholds and epsilon combinations are all invalid")

      epsilon = eugrid[which.min(nllheu), 1]
      ul = eugrid[which.min(nllheu), 2]
      ur = eugrid[which.min(nllheu), 3]

    } else {
      if ((ulseq + epsilon) >= (urseq - epsilon)) stop("lower transition range cannot overlap the upper transition range")
      epsilon = eseq
      ul = ulseq
      ur = urseq
    }
    
    if (fixedeu) { # threshold fixed
      if (is.null(pvector)) {
        pvector[1] = mean(x, trim = 0.2)
        pvector[2] = sd(x)
        initfgpd = fgpd(-x, -ul, std.err = FALSE)
        pvector[3] = initfgpd$sigmau
        pvector[4] = initfgpd$xi
        initfgpd = fgpd(x, ur, std.err = FALSE)
        pvector[5] = initfgpd$sigmau
        pvector[6] = initfgpd$xi
      }
    } else { # threshold as initial value in usual MLE
      if (is.null(pvector)) {
        pvector[1] = mean(x, trim = 0.2)
        pvector[2] = sd(x)
        pvector[3] = epsilon
        pvector[4] = ul
        initfgpd = fgpd(-x, -pvector[4], std.err = FALSE)
        pvector[5] = initfgpd$sigmau
        pvector[6] = initfgpd$xi
        pvector[7] = ur
        initfgpd = fgpd(x, pvector[7], std.err = FALSE)
        pvector[8] = initfgpd$sigmau
        pvector[9] = initfgpd$xi
      } else {
        pvector[9] = pvector[6] # shift upper tail GPD scale and shape to add in ur
        pvector[8] = pvector[5]
        pvector[7] = ur
        pvector[6] = pvector[4] # shift lower tail GPD scale and shape to add in ul
        pvector[5] = pvector[3]
        pvector[4] = ul
        pvector[3] = epsilon
      }
    }
  }
  if (fixedeu) { # fixed threshold (separable) likelihood
    nllh = nleuitmgng(pvector, epsilon, ul, ur, x)
    
    if (is.infinite(nllh)) {
      pvector[4] = 0.1
      pvector[6] = 0.1
      nllh = nleuitmgng(pvector, epsilon, ul, ur, x)
    }
    
    if (is.infinite(nllh)) stop("initial parameter values are invalid")
  
    fit = optim(par = as.vector(pvector), fn = nleuitmgng, epsilon = epsilon, ul = ul, ur = ur, x = x,
      finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)    
    
    nmean = fit$par[1]
    nsd = fit$par[2]
    sigmaul = fit$par[3]
    xil = fit$par[4]
    sigmaur = fit$par[5]
    xir = fit$par[6]
    
  } else { # complete (non-separable) likelihood
    
    nllh = nlitmgng(pvector, x)
    
    if (is.infinite(nllh)) {
      pvector[6] = 0.1
      pvector[9] = 0.1
      nllh = nlitmgng(pvector, x)
    }
    
    if (is.infinite(nllh)) stop("initial parameter values are invalid")
  
    fit = optim(par = as.vector(pvector), fn = nlitmgng, x = x,
      finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)    
    
    nmean = fit$par[1]
    nsd = fit$par[2]
    epsilon = fit$par[3]
    ul = fit$par[4]
    sigmaul = fit$par[5]
    xil = fit$par[6]
    ur = fit$par[7]
    sigmaur = fit$par[8]
    xir = fit$par[9]
  }

  conv = TRUE
  if ((fit$convergence != 0) | any(fit$par == pvector) | (abs(fit$value) >= 1e6)) {
    conv = FALSE
    warning("check convergence")
  }
  
  kappa = 1/(2 + pnorm(ur, nmean, nsd) - pnorm(ul, nmean, nsd))

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
  
  if (!exists("nllheu")) nllheu = NULL

  list(call = call, x = as.vector(x),
    init = as.vector(pvector), fixedeu = fixedeu, ulseq = ulseq, urseq = urseq, eseq = eseq, nllheuseq = nllheu,
    optim = fit, conv = conv, cov = invhess, mle = fit$par, se = se,
    nllh = fit$value, n = n, nmean = nmean, nsd = nsd, epsilon = epsilon,
    ul = ul, sigmaul = sigmaul, xil = xil, ur = ur, sigmaur = sigmaur, xir = xir, kappa = kappa)
}

#' @export
#' @aliases fitmgng litmgng nlitmgng profluitmgng nluitmgng
#' @rdname  fitmgng

# log-likelihood function for normal bulk with GPD for both tails
# interval transition mixture model
litmgng <- function(x, nmean = 0, nsd = 1, epsilon = nsd,
  ul = 0, sigmaul = 1, xil = 0, ur = 0, sigmaur = 1, xir = 0, log = TRUE) {

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(nmean)
  check.param(nsd)
  check.param(epsilon)
  check.param(ul)                       
  check.param(sigmaul)
  check.param(xil)
  check.param(ur)                       
  check.param(sigmaur)
  check.param(xir)
  check.logic(log)
  
  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    x = x[is.finite(x)] # ignore missing and infinite cases
  }

  check.quant(x)
  n = length(x)
  
  check.inputn(c(length(nmean), length(nsd), length(epsilon), 
    length(ul), length(sigmaul), length(xil), length(ur), length(sigmaur), length(xir)), allowscalar = TRUE)

  # assume NA or NaN are irrelevant as entire lower tail is now modelled
  # inconsistent with evd library definition
  # hence use which() to ignore these

  xul = x[which(x < (ul - epsilon))]
  xl = x[which((x >= (ul - epsilon)) & (x <= (ul + epsilon)))]
  xb = x[which((x > (ul + epsilon)) & (x < (ur - epsilon)))]
  xr = x[which((x >= (ur - epsilon)) & (x <= (ur + epsilon)))]
  xur = x[which(x > (ur + epsilon))]
  nul = length(xul)
  nl = length(xl)
  nb = length(xb)
  nr = length(xr)
  nur = length(xur)

  xit = c(xl, xr)
  nit = nl + nr
  
  if ((nsd <= 0) | (epsilon < 0) | (sigmaul <= 0) | (sigmaur <= 0) |
        ((ul - epsilon) <= min(x)) | ((ur + epsilon) >= max(x)) |
        ((ul + epsilon) >= (ur - epsilon))) {
    l = -Inf
  } else {
    kappa = 1/(2 + pnorm(ur, nmean, nsd) - pnorm(ul, nmean, nsd))
    
    if (n != (nul + nl + nb + nr + nur)) {
      stop("total non-finite sample size is not equal to those above/below/between intervals or within it")    
    }

    syul = 1 + xil * (ul - xul) / sigmaul
    syur = 1 + xir * (xur - ur) / sigmaur  
    yb = (xb - nmean) / nsd    # used for normal
  
    if ((min(syul) <= 0) | (min(syur) <= 0)) {
      l = -Inf
    } else { 
       l = lgpd(xur, ur, sigmaur, xir)
       l = l + lgpd(-xul, -ul, sigmaul, xil)
       l = l - nb * log(2 * pi * nsd ^ 2) / 2 - sum(yb ^ 2) / 2
       if (nit > 0) l = l + sum(ditmgng(xit, nmean, nsd, epsilon, ul, sigmaul, xil,
                                                                      ur, sigmaur, xir, log = TRUE))
       l = l + (n - nit) * log(kappa)
    }
  }

  if (!log) l = exp(l)
  
  l
}

#' @export
#' @aliases fitmgng litmgng nlitmgng profluitmgng nluitmgng
#' @rdname  fitmgng

# negative log-likelihood function for normal bulk with GPD for both tails
# interval transition mixture model
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlitmgng <- function(pvector, x, finitelik = FALSE) {

  np = 9 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.logic(finitelik)

  nmean = pvector[1]
  nsd = pvector[2]
  epsilon = pvector[3]
  ul = pvector[4]
  sigmaul = pvector[5]
  xil = pvector[6]
  ur = pvector[7]
  sigmaur = pvector[8]
  xir = pvector[9]

  nllh = -litmgng(x, nmean, nsd, epsilon, ul, sigmaul, xil, ur, sigmaur, xir) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}

#' @export
#' @aliases fitmgng litmgng nlitmgng profluitmgng nluitmgng
#' @rdname  fitmgng

# profile negative log-likelihood function for given threshold for
# normal bulk with GPD for both tails interval transition mixture model
# designed for apply to loop over vector of epsilons and thresholds (hence c(ul, ur) vector is first input)
profleuitmgng <- function(eulr, pvector, x,
  method = "BFGS", control = list(maxit = 10000), finitelik = TRUE, ...) {

  np = 9 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np - 3, allownull = TRUE)
  check.param(eulr, allowvec = TRUE)
  check.nparam(eulr, nparam = 3)
  check.posparam(eulr[1], allowzero = TRUE)
  check.param(eulr[2])
  check.param(eulr[3])
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.optim(method)
  check.control(control)
  check.logic(finitelik)

  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    x = x[is.finite(x)] # ignore missing and infinite cases
  }

  check.quant(x)
  
  # check initial values for other parameters, try usual alternative
  if (!is.null(pvector)) {
    nllh = nleuitmgng(pvector, eulr[1], eulr[2], eulr[3], x)
    
    if (is.infinite(nllh)) pvector = NULL
  }

  if (is.null(pvector)) {
    pvector[1] = mean(x, trim = 0.2)
    pvector[2] = sd(x)
    initfgpd = fgpd(-x, -eulr[2], std.err = FALSE)
    pvector[3] = initfgpd$sigmau
    pvector[4] = initfgpd$xi
    initfgpd = fgpd(x, eulr[3], std.err = FALSE)
    pvector[5] = initfgpd$sigmau
    pvector[6] = initfgpd$xi
    nllh = nleuitmgng(pvector, eulr[1], eulr[2], eulr[3], x)
  }

  if (is.infinite(nllh)) {
    pvector[4] = 0.1
    pvector[6] = 0.1
    nllh = nleuitmgng(pvector, eulr[1], eulr[2], eulr[3], x)
  }

  # if still invalid then output cleanly
  if (is.infinite(nllh)) {
    warning(paste("initial parameter values for thresholds ul =", eulr[2], "and ur =", eulr[3],
                  "and epsilon=", eulr[1], "are invalid"))
    fit = list(par = rep(NA, np), value = Inf, counts = 0, convergence = NA, 
      message = "initial values invalid", hessian = rep(NA, np))
  } else {

    fit = optim(par = as.vector(pvector), fn = nleuitmgng, epsilon = eulr[1], ul = eulr[2], ur = eulr[3], x = x,
      finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)
  }
    
  if (finitelik & is.infinite(fit$value)) {
    fit$value = sign(fit$value) * 1e6
  }

  fit$value
}

#' @export
#' @aliases fitmgng litmgng nlitmgng profluitmgng nluitmgng
#' @rdname  fitmgng

# negative log-likelihood function for normal bulk with GPD for both tails
# interval transition mixture model
# (wrapper for likelihood, designed for threshold to be fixed and other parameters optimised)
nleuitmgng <- function(pvector, epsilon, ul, ur, x, finitelik = FALSE) {

  np = 9 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np - 3)
  check.posparam(epsilon, allowzero = TRUE)
  check.param(ul)
  check.param(ur)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.logic(finitelik)
    
  nmean = pvector[1]
  nsd = pvector[2]
  sigmaul = pvector[3]
  xil = pvector[4]
  sigmaur = pvector[5]
  xir = pvector[6]

  nllh = -litmgng(x, nmean, nsd, epsilon, ul, sigmaul, xil, ur, sigmaur, xir)
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}
