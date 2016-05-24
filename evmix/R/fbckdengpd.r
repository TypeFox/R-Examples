#' @export
#' 
#' @title MLE Fitting of Boundary Corrected Kernel Density Estimate for Bulk and GPD Tail Extreme Value Mixture Model
#'
#' @description Maximum likelihood estimation for fitting the extreme value 
#' mixture model with boundary corrected kernel density estimate for bulk distribution upto the threshold and conditional
#' GPD above threshold. With options for profile likelihood estimation for threshold and
#' fixed threshold approach.
#'
#' @inheritParams fkden
#' @inheritParams fnormgpd
#' @inheritParams dbckdengpd
#' @inheritParams fgpd
#' 
#' @details The extreme value mixture model with boundary corrected kernel density estimate (BCKDE) for bulk and GPD tail is 
#' fitted to the entire dataset using maximum likelihood estimation. The estimated
#' parameters, variance-covariance matrix and their standard errors are automatically
#' output.
#' 
#' See help for \code{\link[evmix:fnormgpd]{fnormgpd}} for details, type \code{help fnormgpd}. 
#' Only the different features are outlined below for brevity.
#' 
#' The full parameter vector is
#' (\code{lambda}, \code{u}, \code{sigmau}, \code{xi}) if threshold is also estimated and
#' (\code{lambda}, \code{sigmau}, \code{xi}) for profile likelihood or fixed threshold approach.
#' 
#' Negative data are ignored.
#' 
#' Cross-validation likelihood is used for BCKDE, but standard likelihood is used
#' for GPD component. See help for \code{\link[evmix:fkden]{fkden}} for details,
#' type \code{help fkden}.
#' 
#' The alternate bandwidth definitions are discussed in the 
#' \code{\link[evmix:kernels]{kernels}}, with the \code{lambda} as the default
#' used in the likelihood fitting. The \code{bw} specification is the same as
#' used in the \code{\link[stats:density]{density}} function.
#' 
#' The possible kernels are also defined in \code{\link[evmix:kernels]{kernels}}
#' with the \code{"gaussian"} as the default choice.
#' 
#' Unlike the standard KDE, there is no general rule-of-thumb bandwidth for all these
#' estimators, with only certain methods having a guideline in the literature, so none
#' have been implemented. Hence, a bandwidth must always be specified.
#' 
#' The \code{simple}, \code{renorm}, \code{beta1}, \code{beta2} \code{gamma1} and \code{gamma2}
#' boundary corrected kernel density estimates require renormalisation, achieved
#' by numerical integration, so are very time consuming.
#' 
#' @section Boundary Correction Methods:
#' 
#' See \code{\link[evmix:bckden]{dbckden}} for details of BCKDE methods.
#' 
#' @section Warning:
#' See important warnings about cross-validation likelihood estimation in 
#' \code{\link[evmix:fkden]{fkden}}, type \code{help fkden}.
#' 
#' See important warnings about boundary correction approaches in 
#' \code{\link[evmix:bckden]{dbckden}}, type \code{help bckden}.
#' 
#' @return \code{\link[evmix:fbckdengpd]{lbckdengpd}}, \code{\link[evmix:fbckdengpd]{nlbckdengpd}},
#' and \code{\link[evmix:fbckdengpd]{nlubckdengpd}} give the log-likelihood,
#' negative log-likelihood and profile likelihood for threshold. Profile likelihood
#' for single threshold is given by \code{\link[evmix:fbckdengpd]{proflubckdengpd}}.
#' \code{\link[evmix:fbckdengpd]{fbckdengpd}} returns a simple list with the following elements
#'
#' \tabular{ll}{
#'  \code{call}:      \tab \code{optim} call\cr
#'  \code{x}:         \tab data vector \code{x}\cr
#'  \code{init}:      \tab \code{pvector}\cr
#'  \code{fixedu}:    \tab fixed threshold, logical\cr
#'  \code{useq}:      \tab threshold vector for profile likelihood or scalar for fixed threshold\cr
#'  \code{nllhuseq}:  \tab profile negative log-likelihood at each threshold in useq\cr
#'  \code{optim}:     \tab complete \code{optim} output\cr
#'  \code{mle}:       \tab vector of MLE of parameters\cr
#'  \code{cov}:       \tab variance-covariance matrix of MLE of parameters\cr
#'  \code{se}:        \tab vector of standard errors of MLE of parameters\cr
#'  \code{rate}:      \tab \code{phiu} to be consistent with \code{\link[evd:fpot]{evd}}\cr
#'  \code{nllh}:      \tab minimum negative log-likelihood\cr
#'  \code{n}:         \tab total sample size\cr
#'  \code{lambda}:    \tab MLE of lambda (kernel half-width)\cr
#'  \code{u}:         \tab threshold (fixed or MLE)\cr
#'  \code{sigmau}:    \tab MLE of GPD scale\cr
#'  \code{xi}:        \tab MLE of GPD shape\cr
#'  \code{phiu}:      \tab MLE of tail fraction (bulk model or parameterised approach)\cr
#'  \code{se.phiu}:   \tab standard error of MLE of tail fraction\cr
#'  \code{bw}:        \tab MLE of bw (kernel standard deviations)\cr
#'  \code{kernel}:    \tab kernel name\cr
#'  \code{bcmethod}:  \tab boundary correction method\cr
#'  \code{proper}:    \tab logical, whether renormalisation is requested\cr
#'  \code{nn}:        \tab non-negative correction method\cr
#'  \code{offset}:    \tab offset for log transformation method\cr
#'  \code{xmax}:      \tab maximum value of scaled beta or copula
#' }
#' 
#' @note 
#' See notes in \code{\link[evmix:fnormgpd]{fnormgpd}} for details, type \code{help fnormgpd}.
#' Only the different features are outlined below for brevity.
#' 
#' No default initial values for parameter vector are provided, so will stop evaluation if
#' \code{pvector} is left as \code{NULL}. Avoid setting the starting value for the shape parameter to
#' \code{xi=0} as depending on the optimisation method it may be get stuck.
#' 
#' The data and kernel centres are both vectors. Infinite, missing and negative sample values
#' (and kernel centres) are dropped.
#' 
#' @references
#' \url{http://www.math.canterbury.ac.nz/~c.scarrott/evmix}
#' 
#' \url{http://en.wikipedia.org/wiki/Kernel_density_estimation}
#' 
#' \url{http://en.wikipedia.org/wiki/Cross-validation_(statistics)}
#' 
#' \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
#' 
#' Scarrott, C.J. and MacDonald, A. (2012). A review of extreme value
#' threshold estimation and uncertainty quantification. REVSTAT - Statistical
#' Journal 10(1), 33-59. Available from \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}
#' 
#' Hu, Y. (2013). Extreme value mixture modelling: An R package and simulation study.
#' MSc (Hons) thesis, University of Canterbury, New Zealand.
#' \url{http://ir.canterbury.ac.nz/simple-search?query=extreme&submit=Go}
#' 
#' Bowman, A.W. (1984). An alternative method of cross-validation for the smoothing of
#' density estimates. Biometrika 71(2), 353-360.
#' 
#' Duin, R.P.W. (1976). On the choice of smoothing parameters for Parzen estimators of
#' probability density functions. IEEE Transactions on Computers C25(11), 1175-1179.
#' 
#' MacDonald, A., Scarrott, C.J., Lee, D., Darlow, B., Reale, M. and Russell, G. (2011).
#' A flexible extreme value mixture model. Computational Statistics and Data Analysis
#' 55(6), 2137-2157.
#' 
#' MacDonald, A., C. J. Scarrott, and D. S. Lee (2011). Boundary correction, consistency
#' and robustness of kernel densities using extreme value theory. Submitted.
#' Available from: \url{http://www.math.canterbury.ac.nz/~c.scarrott}.
#' 
#' Wand, M. and Jones, M.C. (1995). Kernel Smoothing. Chapman && Hall.
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @section Acknowledgments: See Acknowledgments in
#'   \code{\link[evmix:fnormgpd]{fnormgpd}}, type \code{help fnormgpd}. Based on code
#' by Anna MacDonald produced for MATLAB.
#' 
#' @seealso \code{\link[evmix:kernels]{kernels}}, \code{\link[evmix:kfun]{kfun}},
#'  \code{\link[stats:density]{density}}, \code{\link[stats:bandwidth]{bw.nrd0}}
#' and \code{\link[ks:kde.1d]{dkde}} in \code{\link[ks:kde.1d]{ks}} package.
#'  \code{\link[evmix:fgpd]{fgpd}} and \code{\link[evmix:gpd]{gpd}}
#' @aliases fbckdengpd lbckdengpd nlbckdengpd proflubckdengpd nlubckdengpd
#' @family  kdengpd kdengpdcon fkdengpd fkdengpdcon normgpd fnormgpd
#'          kden bckden bckdengpd bckdengpdcon fkden fbckden fbckdengpd fbckdengpdcon
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' par(mfrow = c(2, 1))
#' 
#' x = rgamma(500, 2, 1)
#' xx = seq(-0.1, 10, 0.01)
#' y = dgamma(xx, 2, 1)
#' 
#' # Bulk model based tail fraction
#' pinit = c(0.1, quantile(x, 0.9), 1, 0.1) # initial values required for BCKDE
#' fit = fbckdengpd(x, pvector = pinit, bcmethod = "cutnorm")
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-0.1, 10))
#' lines(xx, y)
#' with(fit, lines(xx, dbckdengpd(xx, x, lambda, u, sigmau, xi, bcmethod = "cutnorm"), col="red"))
#' abline(v = fit$u, col = "red")
#'   
#' # Parameterised tail fraction
#' fit2 = fbckdengpd(x, phiu = FALSE, pvector = pinit, bcmethod = "cutnorm")
#' with(fit2, lines(xx, dbckdengpd(xx, x, lambda, u, sigmau, xi, phiu, bc = "cutnorm"), col="blue"))
#' abline(v = fit2$u, col = "blue")
#' legend("topright", c("True Density","Bulk Tail Fraction","Parameterised Tail Fraction"),
#'   col=c("black", "red", "blue"), lty = 1)
#'   
#' # Profile likelihood for initial value of threshold and fixed threshold approach
#' pinit = c(0.1, 1, 0.1) # notice threshold dropped from initial values
#' fitu = fbckdengpd(x, useq = seq(1, 6, length = 20), pvector = pinit, bcmethod = "cutnorm")
#' fitfix = fbckdengpd(x, useq = seq(1, 6, length = 20), fixedu = TRUE, pv = pinit, bc = "cutnorm")
#' 
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-0.1, 10))
#' lines(xx, y)
#' with(fit, lines(xx, dbckdengpd(xx, x, lambda, u, sigmau, xi, bc = "cutnorm"), col="red"))
#' abline(v = fit$u, col = "red")
#' with(fitu, lines(xx, dbckdengpd(xx, x, lambda, u, sigmau, xi, bc = "cutnorm"), col="purple"))
#' abline(v = fitu$u, col = "purple")
#' with(fitfix, lines(xx, dbckdengpd(xx, x, lambda, u, sigmau, xi, bc = "cutnorm"), col="darkgreen"))
#' abline(v = fitfix$u, col = "darkgreen")
#' legend("topright", c("True Density","Default initial value (90% quantile)",
#'  "Prof. lik. for initial value", "Prof. lik. for fixed threshold"),
#'  col=c("black", "red", "purple", "darkgreen"), lty = 1)
#' }
#'

# maximum likelihood fitting for boundary corrected kernel density estimate for bulk with GPD for upper tail
fbckdengpd <- function(x, phiu = TRUE, useq = NULL, fixedu = FALSE, pvector = NULL, kernel = "gaussian",
  bcmethod = "simple", proper = TRUE, nn = "jf96", offset = NULL, xmax = NULL,
  add.jitter = FALSE, factor = 0.1, amount = NULL,
  std.err = TRUE, method = "BFGS", control = list(maxit = 10000), finitelik = TRUE, ...) {

  call <- match.call()
    
  np = 4 # maximum number of parameters

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.logic(phiu)
  check.posparam(useq, allowvec = TRUE, allownull = TRUE)
  check.logic(fixedu)
  check.logic(std.err)
  check.optim(method)
  check.control(control)
  check.logic(finitelik)

  check.kernel(kernel)
  check.bcmethod(bcmethod)
  check.logic(proper)
  check.nn(nn)
  check.offset(offset, bcmethod, allowzero = TRUE)
  check.posparam(xmax, allownull = TRUE)  
  check.posparam(factor)
  check.posparam(amount, allownull = TRUE)
  check.logic(add.jitter)
  
  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    x = x[is.finite(x)] # ignore missing and infinite cases
  }

  if (any(x < 0)) {
    warning("negative values have been removed")
    x = x[x >= 0]
  }
  
  check.quant(x)
  n = length(x)

  if (add.jitter) x = pmax(jitter(x, factor, amount), 0)

  xuniq = unique(x)
  if (length(xuniq) < (0.95*n))
    warning("data may be rounded, as more than 5% are ties, so bandwidth could be biased to zero")

  if (is.null(pvector)) {
    stop("Initial values for parameter vector must be provided")
  }

  linit = pvector[1]
  check.posparam(linit)
  
  if ((bcmethod == "copula") & (linit >= 1))
    stop("bandwidth must between (0, 1) for copula method")  
    
  upboundmethods = c("beta1", "beta2", "copula")
  if (!is.null(xmax) & !(bcmethod %in% upboundmethods))
    warning(paste("xmax only relevant for boundary correction methods", upboundmethods, collapse = " "))
  
  if (bcmethod %in% upboundmethods) {
    if (is.null(xmax)) stop("xmax is NULL")
    
    if (max(x) > xmax) stop("largest kernel centre must be below xmax")

    if (any(x == 0)) {
      warning("kernel centres of zero are invalid for gamma or beta method so ignored")
      x = x[x != 0]
    }

    if ((bcmethod != "gamma1") & (bcmethod != "gamma2")) {
      if (any(x == xmax)) {
        warning("kernel centres of xmax are invalid for beta or copula method so ignored")
        x = x[x != xmax]
      }
    }
    # need to recheck there are some valid kernel centres after these exclusions
    check.quant(x)
    n = length(x)

    if (max(x) > xmax) stop("largest kernel centre must be below xmax")   
  }

  # It is not always easy to choose a sensible good initial value for lambda
  # So try adjusting lambda up and down a little to find a valid one to start off
  llhinit = lbckden(x, lambda = linit, kernel = kernel,
    bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax)

  tryi = 0
  lfirst = linit
  
  # try upto 2^5 larger than original
  while (is.infinite(llhinit) & (tryi < 5)) {
    linit = linit*2
    llhinit = lbckden(x, lambda = linit, kernel = kernel,
      bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax)
    tryi = tryi + 1
  }
  
  # try upto 2^-5 smaller than original
  if (is.infinite(llhinit)) {
    tryi = 0
    linit = lfirst
    while (is.infinite(llhinit) & (tryi < 5)) {
      linit = linit/2
      llhinit = lbckden(x, lambda = linit, kernel = kernel, 
        bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax)
      tryi = tryi + 1
    }
  }
  
  if (is.infinite(llhinit))
    stop("likelihood is undefined for initial bandwidth try another value")  

  if (tryi != 0)
    warning(paste("initial bandwidth was invalid, so linit=", linit, "is used instead"))
  
  # add back the initial value for bandwidth
  pvector[1] = linit
  
  if ((method == "L-BFGS-B") | (method == "BFGS")) finitelik = TRUE
  
  # useq must be specified if threshold is fixed
  if (fixedu & is.null(useq))
    stop("for fixed threshold approach, useq must be specified (as scalar or vector)")
  
  # Check if profile likelihood or fixed threshold is being used
  # and determine initial values for parameters in each case
  if (is.null(useq)) { # not profile or fixed

    check.nparam(pvector, nparam = np, allownull = TRUE)
        
  } else { # profile or fixed
    
    check.nparam(pvector, nparam = np - 1, allownull = TRUE)

    # profile likelihood for threshold or scalar given
    if (length(useq) != 1) {
      
      # remove thresholds with less than 5 excesses
      useq = useq[sapply(useq, FUN = function(u, x) sum(x > u) > 5, x = x)]
      check.posparam(useq, allowvec = TRUE)
      
      nllhu = sapply(useq, proflubckdengpd, pvector = pvector, x = x, phiu = phiu, kernel = kernel,
        bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax,
        method = method, control = control, finitelik = finitelik, ...)
      
      if (all(!is.finite(nllhu))) stop("thresholds are all invalid")
      u = useq[which.min(nllhu)]

    } else {
      u = useq
    }

    if (!fixedu) { # threshold as initial value in usual MLE
      pvector[4] = pvector[3] # shift GPD scale and shape to add in u
      pvector[3] = pvector[2]
      pvector[2] = u
    }
  }

  if (fixedu) { # fixed threshold (separable) likelihood
    nllh = nlubckdengpd(pvector, u, x, phiu, kernel = kernel,
                        bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax)
    if (is.infinite(nllh)) {
      pvector[3] = 0.1
      nllh = nlubckdengpd(pvector, u, x, phiu, kernel = kernel,
                          bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax)
    }
    if (is.infinite(nllh)) stop("initial parameter values are invalid")

    fit = optim(par = as.vector(pvector), fn = nlubckdengpd, u = u, x = x, phiu = phiu, kernel = kernel,
      bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax,
      finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)    
    
    lambda = fit$par[1]
    sigmau = fit$par[2]
    xi = fit$par[3]
    
  } else { # complete (non-separable) likelihood
    
    nllh = nlbckdengpd(pvector, x, phiu, kernel = kernel,
                       bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax)
    if (is.infinite(nllh)) {
      pvector[4] = 0.1
      nllh = nlbckdengpd(pvector, x, phiu, kernel = kernel,
                         bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax)
    }
    if (is.infinite(nllh)) stop("initial parameter values are invalid")
  
    fit = optim(par = as.vector(pvector), fn = nlbckdengpd, x = x, phiu = phiu, kernel = kernel,
      bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax,
      finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)    
    
    lambda = fit$par[1]
    u = fit$par[2]
    sigmau = fit$par[3]
    xi = fit$par[4]
  }

  kernelmethods = c("simple", "cutnorm", "renorm", "reflect", "logtrans")
  if (bcmethod %in% kernelmethods) {
    bw = kbw(lambda, kernel)
  } else {
    bw = NA
  }

  conv = TRUE
  if ((fit$convergence != 0) | any(fit$par == pvector) | (abs(fit$value) >= 1e6)) {
    conv = FALSE
    warning("check convergence")
  }

  pu = pbckden(u, x, lambda, kernel = kernel,
    bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax)
  if (phiu) {
    phiu = 1 - pu
    se.phiu = NA
  } else {
    phiu = mean(x > u, na.rm = TRUE)
    se.phiu = sqrt(phiu * (1 - phiu) / n)
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

  list(call = call, x = as.vector(x),
    init = as.vector(pvector), fixedu = fixedu, useq = useq, nllhuseq = nllhu,
    optim = fit, conv = conv, cov = invhess, mle = fit$par, se = se, rate = phiu,
    nllh = fit$value, n = n,
    lambda = lambda, u = u, sigmau = sigmau, xi = xi, phiu = phiu, se.phiu = se.phiu, bw = bw, kernel = kernel,
    bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax)
}

#' @export
#' @aliases fbckdengpd lbckdengpd nlbckdengpd proflubckdengpd nlubckdengpd
#' @rdname  fbckdengpd

# log-likelihood function for boundary corrected kernel density estimate for bulk with GPD for upper tail
# cross-validation for KDE component
lbckdengpd <- function(x, lambda = NULL, u = 0, sigmau = 1, xi = 0, phiu = TRUE,
  bw = NULL, kernel = "gaussian",
  bcmethod = "simple", proper = TRUE, nn = "jf96", offset = NULL, xmax = NULL, log = TRUE) {

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(lambda, allownull = TRUE)
  check.param(bw, allownull = TRUE)
  check.param(u)
  check.param(sigmau)
  check.param(xi)
  check.phiu(phiu, allowfalse = TRUE)
  check.logic(log)

  check.kernel(kernel)
  check.bcmethod(bcmethod)
  check.logic(proper)
  check.nn(nn)
  check.offset(offset, bcmethod, allowzero = TRUE)
  check.posparam(xmax, allownull = TRUE)  

  kernel = ifelse(kernel == "rectangular", "uniform", kernel)
  kernel = ifelse(kernel == "normal", "gaussian", kernel)
  
  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    x = x[is.finite(x)] # ignore missing and infinite cases
  }

  if (any(x < 0)) {
    warning("negative values have been removed")
    x = x[x >= 0]
  }

  check.quant(x)
  n = length(x)

  # if bcmethod does not use standard kernels then lambda must be specified
  # then bw can be used, but lambda should be defaulted to if available
  kernelmethods = c("simple", "cutnorm", "renorm", "reflect", "logtrans")
  if (!(bcmethod %in% kernelmethods)) {
    if (is.null(lambda))
      stop(paste("bandwidth bw only relevant for", kernelmethods, collapse = " "))
  } else {
    if (is.null(lambda) & is.null(bw)) stop("lambda and bw cannot both be NULL")
    
    if (is.null(lambda)) lambda = klambda(bw, kernel)
  }
  
  check.inputn(c(length(lambda), length(u), length(sigmau), length(xi), length(phiu)), allowscalar = TRUE)

  upboundmethods = c("beta1", "beta2", "copula")
  if (!is.null(xmax) & !(bcmethod %in% upboundmethods))
    warning(paste("xmax only relevant for boundary correction methods", upboundmethods, collapse = " "))
  
  if (bcmethod %in% upboundmethods) {
    if (is.null(xmax)) stop("xmax is NULL")
    
    if (max(x) > xmax) stop("largest kernel centre must be below xmax")

    if (any(x == 0)) {
      warning("kernel centres of zero are invalid for gamma or beta method so ignored")
      x = x[x != 0]
    }

    if ((bcmethod != "gamma1") & (bcmethod != "gamma2")) {
      if (any(x == xmax)) {
        warning("kernel centres of xmax are invalid for beta or copula method so ignored")
        x = x[x != xmax]
      }
    }
    # need to recheck there are some valid kernel centres after these exclusions
    check.quant(x)
    n = length(x)
  }

  # assume NA or NaN are irrelevant as entire lower tail is now modelled
  # inconsistent with evd library definition
  # hence use which() to ignore these

  xu = x[which(x > u)]
  nu = length(xu)
  xb = x[which(x <= u)]
  nb = length(xb)

  if (n != nb + nu) {
    stop("total non-finite sample size is not equal to those above threshold and those below or equal to it")    
  }

  if ((lambda <= 0) | ((bcmethod == "copula") & (lambda >= 1)) |
      ((bcmethod == "beta1") & (lambda >= 0.25*ifelse(is.null(xmax), Inf, xmax))) | 
      ((bcmethod == "beta2") & (lambda >= 0.25*ifelse(is.null(xmax), Inf, xmax))) |
      (sigmau <= 0) | (u <= min(x)) | (u >= max(x))) {
    l = -Inf
  } else {
    if (is.logical(phiu)) {
      pu = pbckden(u, x, lambda, kernel = kernel, bcmethod = bcmethod,
        proper = proper, nn = nn, offset = offset, xmax = xmax)
      if (phiu) {
        phiu = 1 - pu
      } else {
        phiu = nu / n
      }
    }
    phib = (1 - phiu) / pu
  
    syu = 1 + xi * (xu - u) / sigmau  
  
    if ((min(syu) <= 0) | (phiu <= 0) | (phiu >= 1) | (pu <= 0) | (pu >= 1) | ifelse(is.null(xmax), FALSE, u >= xmax)) {
      l = -Inf
    } else { 
      l = lgpd(xu, u, sigmau, xi, phiu)
      l = l + lbckden(xb, lambda, kernel = kernel, extracentres = xu, bcmethod = bcmethod,
        proper = proper, nn = nn, offset = offset, xmax = xmax, log = TRUE) + nb*log(phib)
    }
  }
  
  if (!log) l = exp(l)
  
  l
}

#' @export
#' @aliases fbckdengpd lbckdengpd nlbckdengpd proflubckdengpd nlubckdengpd
#' @rdname  fbckdengpd

# negative log-likelihood function for boundary corrected kernel density estimate for bulk with GPD for upper tail
# cross-validation for KDE component
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlbckdengpd <- function(pvector, x, phiu = TRUE, kernel = "gaussian",
  bcmethod = "simple", proper = TRUE, nn = "jf96", offset = NULL, xmax = NULL, finitelik = FALSE) {

  np = 4 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.phiu(phiu, allowfalse = TRUE)
  check.logic(finitelik)

  check.kernel(kernel)
  check.bcmethod(bcmethod)
  check.logic(proper)
  check.nn(nn)
  check.offset(offset, bcmethod, allowzero = TRUE)
  check.posparam(xmax, allownull = TRUE)  

  kernel = ifelse(kernel == "rectangular", "uniform", kernel)
  kernel = ifelse(kernel == "normal", "gaussian", kernel)

  lambda = pvector[1]
  u = pvector[2]
  sigmau = pvector[3]
  xi = pvector[4]

  nllh = -lbckdengpd(x, lambda, u, sigmau, xi, phiu, kernel = kernel,
    bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax)
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}

#' @export
#' @aliases fbckdengpd lbckdengpd nlbckdengpd proflubckdengpd nlubckdengpd
#' @rdname  fbckdengpd

# profile negative log-likelihood function for given threshold for
# boundary corrected kernel density estimate for bulk with GPD for upper tail
# designed for sapply to loop over vector of thresholds (hence u is first input)
# cross-validation for KDE component
proflubckdengpd <- function(u, pvector, x, phiu = TRUE, kernel = "gaussian",
  bcmethod = "simple", proper = TRUE, nn = "jf96", offset = NULL, xmax = NULL,
  method = "BFGS", control = list(maxit = 10000), finitelik = TRUE, ...) {

  np = 4 # maximum number of parameters
  
  # Check properties of inputs
  check.nparam(pvector, nparam = np - 1, allownull = TRUE)
  check.posparam(u)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.phiu(phiu, allowfalse = TRUE)
  check.optim(method)
  check.control(control)
  check.logic(finitelik)

  check.kernel(kernel)
  check.bcmethod(bcmethod)
  check.logic(proper)
  check.nn(nn)
  check.offset(offset, bcmethod, allowzero = TRUE)
  check.posparam(xmax, allownull = TRUE)  
  
  kernel = ifelse(kernel == "rectangular", "uniform", kernel)
  kernel = ifelse(kernel == "normal", "gaussian", kernel)

  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    x = x[is.finite(x)] # ignore missing and infinite cases
  }

  if (any(x < 0)) {
    warning("negative values have been removed")
    x = x[x >= 0]
  }

  check.quant(x)
  
  upboundmethods = c("beta1", "beta2", "copula")
  if (!is.null(xmax) & !(bcmethod %in% upboundmethods))
    warning(paste("xmax only relevant for boundary correction methods", upboundmethods, collapse = " "))
  
  if (bcmethod %in% upboundmethods) {
    if (is.null(xmax)) stop("xmax is NULL")
    
    if (max(x) > xmax) stop("largest kernel centre must be below xmax")

    if (any(x == 0)) {
      warning("kernel centres of zero are invalid for gamma or beta method so ignored")
      x = x[x != 0]
    }

    if ((bcmethod != "gamma1") & (bcmethod != "gamma2")) {
      if (any(x == xmax)) {
        warning("kernel centres of xmax are invalid for beta or copula method so ignored")
        x = x[x != xmax]
      }
    }
    # need to recheck there are some valid kernel centres after these exclusions
    check.quant(x)
  }

  if (is.null(pvector)) {
    stop("Initial values for parameter vector must be provided")
  }

  # check initial values for other parameters, try usual alternative
  nllh = nlubckdengpd(pvector, u, x, phiu, kernel = kernel,
                      bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax)
  if (is.infinite(nllh)) {
    pvector[3] = 0.1
    nllh = nlubckdengpd(pvector, u, x, phiu, kernel = kernel,
                        bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax)
  }

  # if still invalid then output cleanly
  if (is.infinite(nllh)) {
    warning(paste("initial parameter values for threshold u =", u, "are invalid"))
    fit = list(par = rep(NA, np), value = Inf, counts = 0, convergence = NA, 
      message = "initial values invalid", hessian = rep(NA, np))
  } else {

    fit = optim(par = as.vector(pvector), fn = nlubckdengpd, u = u, x = x, phiu = phiu, kernel = kernel,
      bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax,
      finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)
  }
    
  if (finitelik & is.infinite(fit$value)) {
    fit$value = sign(fit$value) * 1e6
  }

  fit$value
}

#' @export
#' @aliases fbckdengpd lbckdengpd nlbckdengpd proflubckdengpd nlubckdengpd
#' @rdname  fbckdengpd

# negative log-likelihood function for boundary corrected kernel density estimate for bulk with GPD for upper tail
# (wrapper for likelihood, designed for threshold to be fixed and other parameters optimised)
# cross-validation for KDE component
nlubckdengpd <- function(pvector, u, x, phiu = TRUE, kernel = "gaussian",
  bcmethod = "simple", proper = TRUE, nn = "jf96", offset = NULL, xmax = NULL,
  finitelik = FALSE) {

  np = 4 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np - 1)
  check.posparam(u)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.phiu(phiu, allowfalse = TRUE)
  check.logic(finitelik)

  check.kernel(kernel)
  check.bcmethod(bcmethod)
  check.logic(proper)
  check.nn(nn)
  check.offset(offset, bcmethod, allowzero = TRUE)
  check.posparam(xmax, allownull = TRUE)  
  
  kernel = ifelse(kernel == "rectangular", "uniform", kernel)
  kernel = ifelse(kernel == "normal", "gaussian", kernel)
    
  lambda = pvector[1]
  sigmau = pvector[2]
  xi = pvector[3]

  nllh = -lbckdengpd(x, lambda, u, sigmau, xi, phiu, kernel = kernel,
    bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}
