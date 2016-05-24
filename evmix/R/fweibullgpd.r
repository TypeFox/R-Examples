#' @export
#' 
#' @title MLE Fitting of Weibull Bulk and GPD Tail Extreme Value Mixture Model
#'
#' @description Maximum likelihood estimation for fitting the extreme value 
#' mixture model with Weibull for bulk distribution upto the threshold and conditional
#' GPD above threshold. With options for profile likelihood estimation for threshold and
#' fixed threshold approach.
#'
#' @param wshape  scalar Weibull shape (positive)
#' @param wscale  scalar Weibull scale (positive)
#' @inheritParams fnormgpd
#' @inheritParams fgpd
#' 
#' @details The extreme value mixture model with Weibull bulk and GPD tail is 
#' fitted to the entire dataset using maximum likelihood estimation. The estimated
#' parameters, variance-covariance matrix and their standard errors are automatically
#' output.
#' 
#' See help for \code{\link[evmix:fnormgpd]{fnormgpd}} for details, type \code{help fnormgpd}. 
#' Only the different features are outlined below for brevity.
#' 
#' The full parameter vector is
#' (\code{wshape}, \code{wscale}, \code{u}, \code{sigmau}, \code{xi}) if threshold is also estimated and
#' (\code{wshape}, \code{wscale}, \code{sigmau}, \code{xi}) for profile likelihood or fixed threshold approach.
#' 
#' Non-positive data are ignored (f(0) is infinite for wshape<1).
#' 
#' @return Log-likelihood is given by \code{\link[evmix:fweibullgpd]{lweibullgpd}} and it's
#'   wrappers for negative log-likelihood from \code{\link[evmix:fweibullgpd]{nlweibullgpd}}
#'   and \code{\link[evmix:fweibullgpd]{nluweibullgpd}}. Profile likelihood for single
#'   threshold given by \code{\link[evmix:fweibullgpd]{profluweibullgpd}}. Fitting function
#'   \code{\link[evmix:fweibullgpd]{fweibullgpd}} returns a simple list with the
#'   following elements
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
#'  \code{wshape}:    \tab MLE of Weibull shape\cr
#'  \code{wscale}:    \tab MLE of Weibull scale\cr
#'  \code{u}:         \tab threshold (fixed or MLE)\cr
#'  \code{sigmau}:    \tab MLE of GPD scale\cr
#'  \code{xi}:        \tab MLE of GPD shape\cr
#'  \code{phiu}:      \tab MLE of tail fraction (bulk model or parameterised approach)\cr
#'  \code{se.phiu}:   \tab standard error of MLE of tail fraction\cr
#' }
#' 
#' @note When \code{pvector=NULL} then the initial values are:
#' \itemize{
#'  \item MLE of Weibull parameters assuming entire population is Weibull; and
#'  \item threshold 90\% quantile (not relevant for profile likelihood for threshold or fixed threshold approaches);
#'  \item MLE of GPD parameters above threshold. 
#' }
#' 
#' @references
#' \url{http://www.math.canterbury.ac.nz/~c.scarrott/evmix}
#' 
#' \url{http://en.wikipedia.org/wiki/Weibull_distribution}
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
#' Behrens, C.N., Lopes, H.F. and Gamerman, D. (2004). Bayesian analysis of extreme
#' events with threshold estimation. Statistical Modelling. 4(3), 227-244.
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @section Acknowledgments: See Acknowledgments in
#'   \code{\link[evmix:fnormgpd]{fnormgpd}}, type \code{help fnormgpd}.
#' 
#' @seealso \code{\link[stats:Weibull]{dweibull}},
#'  \code{\link[evmix:fgpd]{fgpd}} and \code{\link[evmix:gpd]{gpd}}
#'  
#' @aliases fweibullgpd lweibullgpd nlweibullgpd profluweibullgpd nluweibullgpd
#' @family  weibullgpd weibullgpdcon fweibullgpd fweibullgpdcon normgpd fnormgpd
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' par(mfrow = c(2, 1))
#' 
#' x = rweibull(1000, shape = 2)
#' xx = seq(-0.1, 4, 0.01)
#' y = dweibull(xx, shape = 2)
#' 
#' # Bulk model based tail fraction
#' fit = fweibullgpd(x)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-0.1, 4))
#' lines(xx, y)
#' with(fit, lines(xx, dweibullgpd(xx, wshape, wscale, u, sigmau, xi), col="red"))
#' abline(v = fit$u, col = "red")
#'   
#' # Parameterised tail fraction
#' fit2 = fweibullgpd(x, phiu = FALSE)
#' with(fit2, lines(xx, dweibullgpd(xx, wshape, wscale, u, sigmau, xi, phiu), col="blue"))
#' abline(v = fit2$u, col = "blue")
#' legend("topright", c("True Density","Bulk Tail Fraction","Parameterised Tail Fraction"),
#'   col=c("black", "red", "blue"), lty = 1)
#'   
#' # Profile likelihood for initial value of threshold and fixed threshold approach
#' fitu = fweibullgpd(x, useq = seq(0.5, 2, length = 20))
#' fitfix = fweibullgpd(x, useq = seq(0.5, 2, length = 20), fixedu = TRUE)
#' 
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-0.1, 4))
#' lines(xx, y)
#' with(fit, lines(xx, dweibullgpd(xx, wshape, wscale, u, sigmau, xi), col="red"))
#' abline(v = fit$u, col = "red")
#' with(fitu, lines(xx, dweibullgpd(xx, wshape, wscale, u, sigmau, xi), col="purple"))
#' abline(v = fitu$u, col = "purple")
#' with(fitfix, lines(xx, dweibullgpd(xx, wshape, wscale, u, sigmau, xi), col="darkgreen"))
#' abline(v = fitfix$u, col = "darkgreen")
#' legend("topright", c("True Density","Default initial value (90% quantile)",
#'  "Prof. lik. for initial value", "Prof. lik. for fixed threshold"),
#'  col=c("black", "red", "purple", "darkgreen"), lty = 1)
#' }
#'   

# maximum likelihood fitting for Weibull bulk with GPD for upper tail
fweibullgpd <- function(x, phiu = TRUE, useq = NULL, fixedu = FALSE, pvector = NULL,
  std.err = TRUE, method = "BFGS", control = list(maxit = 10000), finitelik = TRUE, ...) {

  call <- match.call()
    
  np = 5 # maximum number of parameters

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.logic(phiu)
  check.posparam(useq, allowvec = TRUE, allownull = TRUE)
  check.logic(fixedu)
  check.logic(std.err)
  check.optim(method)
  check.control(control)
  check.logic(finitelik)

  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    x = x[is.finite(x)] # ignore missing and infinite cases
  }

  if (any(x <= 0)) {
    warning("non-positive values have been removed")
    x = x[x > 0]
  }

  check.quant(x)
  n = length(x)

  if ((method == "L-BFGS-B") | (method == "BFGS")) finitelik = TRUE
  
  # useq must be specified if threshold is fixed
  if (fixedu & is.null(useq))
    stop("for fixed threshold approach, useq must be specified (as scalar or vector)")
  
  # Check if profile likelihood or fixed threshold is being used
  # and determine initial values for parameters in each case
  if (is.null(useq)) { # not profile or fixed

    check.nparam(pvector, nparam = np, allownull = TRUE)
    
    if (is.null(pvector)) {
      initfweibull = fitdistr(x, "weibull", lower = c(1e-8, 1e-8))
      pvector[1] = initfweibull$estimate[1]
      pvector[2] = initfweibull$estimate[2]
      pvector[3] = as.vector(quantile(x, 0.9))
      initfgpd = fgpd(x, pvector[3], std.err = FALSE)
      pvector[4] = initfgpd$sigmau
      pvector[5] = initfgpd$xi
    }
    
  } else { # profile or fixed
    
    check.nparam(pvector, nparam = np - 1, allownull = TRUE)

    # profile likelihood for threshold or scalar given
    if (length(useq) != 1) {
      
      # remove thresholds with less than 5 excesses
      useq = useq[sapply(useq, FUN = function(u, x) sum(x > u) > 5, x = x)]
      check.posparam(useq, allowvec = TRUE)
      
      nllhu = sapply(useq, profluweibullgpd, pvector = pvector, x = x, phiu = phiu,
        method = method, control = control, finitelik = finitelik, ...)
      
      if (all(!is.finite(nllhu))) stop("thresholds are all invalid")
      u = useq[which.min(nllhu)]

    } else {
      u = useq
    }

    if (fixedu) { # threshold fixed
      if (is.null(pvector)) {
        initfweibull = fitdistr(x, "weibull", lower = c(1e-8, 1e-8))
        pvector[1] = initfweibull$estimate[1]
        pvector[2] = initfweibull$estimate[2]
        initfgpd = fgpd(x, u, std.err = FALSE)
        pvector[3] = initfgpd$sigmau
        pvector[4] = initfgpd$xi
      }
    } else { # threshold as initial value in usual MLE
      if (is.null(pvector)) {
        initfweibull = fitdistr(x, "weibull", lower = c(1e-8, 1e-8))
        pvector[1] = initfweibull$estimate[1]
        pvector[2] = initfweibull$estimate[2]
        pvector[3] = u
        initfgpd = fgpd(x, pvector[3], std.err = FALSE)
        pvector[4] = initfgpd$sigmau
        pvector[5] = initfgpd$xi
      } else {
        pvector[5] = pvector[4] # shift GPD scale and shape to add in u
        pvector[4] = pvector[3]
        pvector[3] = u
      }
    }
  }

  if (fixedu) { # fixed threshold (separable) likelihood
    nllh = nluweibullgpd(pvector, u, x, phiu)
    if (is.infinite(nllh)) {
      pvector[4] = 0.1
      nllh = nluweibullgpd(pvector, u, x, phiu)    
    }
    if (is.infinite(nllh)) stop("initial parameter values are invalid")
  
    fit = optim(par = as.vector(pvector), fn = nluweibullgpd, u = u, x = x, phiu = phiu,
      finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)    
    
    wshape = fit$par[1]
    wscale = fit$par[2]
    sigmau = fit$par[3]
    xi = fit$par[4]
    
  } else { # complete (non-separable) likelihood
    
    nllh = nlweibullgpd(pvector, x, phiu)
    if (is.infinite(nllh)) {
      pvector[5] = 0.1
      nllh = nlweibullgpd(pvector, x, phiu)    
    }
    if (is.infinite(nllh)) stop("initial parameter values are invalid")
  
    fit = optim(par = as.vector(pvector), fn = nlweibullgpd, x = x, phiu = phiu,
      finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)    
    
    wshape = fit$par[1]
    wscale = fit$par[2]
    u = fit$par[3]
    sigmau = fit$par[4]
    xi = fit$par[5]
  }
  
  conv = TRUE
  if ((fit$convergence != 0) | any(fit$par == pvector) | (abs(fit$value) >= 1e6)) {
    conv = FALSE
    warning("check convergence")
  }

  pu = pweibull(u, wshape, wscale)
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
    wshape = wshape, wscale = wscale, u = u, sigmau = sigmau, xi = xi, phiu = phiu, se.phiu = se.phiu)
}

#' @export
#' @aliases fweibullgpd lweibullgpd nlweibullgpd profluweibullgpd nluweibullgpd
#' @rdname  fweibullgpd

# log-likelihood function for Weibull bulk with GPD for upper tail
lweibullgpd <- function(x, wshape = 1, wscale = 1, u = qweibull(0.9, wshape, wscale),
  sigmau = sqrt(wscale^2 * gamma(1 + 2/wshape) - (wscale * gamma(1 + 1/wshape))^2),
  xi = 0, phiu = TRUE, log = TRUE) {

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(wshape)
  check.param(wscale)
  check.param(u)
  check.param(sigmau)
  check.param(xi)
  check.phiu(phiu, allowfalse = TRUE)
  check.logic(log)

  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    x = x[is.finite(x)] # ignore missing and infinite cases
  }

  if (any(x <= 0)) {
    warning("non-positive values have been removed")
    x = x[x > 0]
  }

  check.quant(x)
  n = length(x)

  check.inputn(c(length(wshape), length(wscale), length(u), length(sigmau), length(xi), length(phiu)),
               allowscalar = TRUE)

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

  if ((wscale <= 0) | (wshape <= 0) | (sigmau <= 0) | (u <= 0) | (u <= min(x)) | (u >= max(x))) {
    l = -Inf
  } else {
    if (is.logical(phiu)) {
      pu = pweibull(u, wshape, wscale)
      if (phiu) {
        phiu = 1 - pu
      } else {
        phiu = nu / n
      }
    }
    phib = (1 - phiu) / pu
  
    syu = 1 + xi * (xu - u) / sigmau  

    if ((min(syu) <= 0) | (phiu <= 0) | (phiu >= 1) | (pu < .Machine$double.eps) | (pu <= 0) | (pu >= 1)) {
      l = -Inf
    } else { 
      l = lgpd(xu, u, sigmau, xi, phiu)
      l = l + (wshape - 1) * sum(log(xb)) - sum(xb^wshape) / wscale^wshape + nb * log(wshape) - nb * wshape * log(wscale)  + nb * log(phib)
    }
  }
  
  if (!log) l = exp(l)
  
  l
}

#' @export
#' @aliases fweibullgpd lweibullgpd nlweibullgpd profluweibullgpd nluweibullgpd
#' @rdname  fweibullgpd

# negative log-likelihood function for Weibull bulk with GPD for upper tail
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlweibullgpd <- function(pvector, x, phiu = TRUE, finitelik = FALSE) {

  np = 5 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.phiu(phiu, allowfalse = TRUE)
  check.logic(finitelik)

  wshape = pvector[1]
  wscale = pvector[2]
  u = pvector[3]
  sigmau = pvector[4]
  xi = pvector[5]

  nllh = -lweibullgpd(x, wshape, wscale, u, sigmau, xi, phiu) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}

#' @export
#' @aliases fweibullgpd lweibullgpd nlweibullgpd profluweibullgpd nluweibullgpd
#' @rdname  fweibullgpd

# profile negative log-likelihood function for given threshold for
# Weibull bulk with GPD for upper tail
# designed for sapply to loop over vector of thresholds (hence u is first input)
profluweibullgpd <- function(u, pvector, x, phiu = TRUE, method = "BFGS",
  control = list(maxit = 10000), finitelik = TRUE, ...) {

  np = 5 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np - 1, allownull = TRUE)
  check.posparam(u)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.phiu(phiu, allowfalse = TRUE)
  check.optim(method)
  check.control(control)
  check.logic(finitelik)

  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    x = x[is.finite(x)] # ignore missing and infinite cases
  }

  if (any(x <= 0)) {
    warning("non-positive values have been removed")
    x = x[x > 0]
  }

  check.quant(x)

  # check initial values for other parameters, try usual alternative
  if (!is.null(pvector)) {
    nllh = nluweibullgpd(pvector, u, x, phiu)
    
    if (is.infinite(nllh)) pvector = NULL
  }

  if (is.null(pvector)) {
    initfweibull = fitdistr(x, "weibull", lower = c(1e-8, 1e-8))
    pvector[1] = initfweibull$estimate[1]
    pvector[2] = initfweibull$estimate[2]
    initfgpd = fgpd(x, u, std.err = FALSE)
    pvector[3] = initfgpd$sigmau
    pvector[4] = initfgpd$xi
    nllh = nluweibullgpd(pvector, u, x, phiu)
  }  

  if (is.infinite(nllh)) {
    pvector[4] = 0.1
    nllh = nluweibullgpd(pvector, u, x, phiu)    
  }

  # if still invalid then output cleanly
  if (is.infinite(nllh)) {
    warning(paste("initial parameter values for threshold u =", u, "are invalid"))
    fit = list(par = rep(NA, np), value = Inf, counts = 0, convergence = NA, 
      message = "initial values invalid", hessian = rep(NA, np))
  } else {

    fit = optim(par = as.vector(pvector), fn = nluweibullgpd, u = u, x = x, phiu = phiu,
    finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)
  }
    
  if (finitelik & is.infinite(fit$value)) {
    fit$value = sign(fit$value) * 1e6
  }

  fit$value
}

#' @export
#' @aliases fweibullgpd lweibullgpd nlweibullgpd profluweibullgpd nluweibullgpd
#' @rdname  fweibullgpd

# negative log-likelihood function for Weibull bulk with GPD for upper tail
# (wrapper for likelihood, designed for threshold to be fixed and other parameters optimised)
nluweibullgpd <- function(pvector, u, x, phiu = TRUE, finitelik = FALSE) {

  np = 5 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np - 1)
  check.posparam(u)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.phiu(phiu, allowfalse = TRUE)
  check.logic(finitelik)
    
  wshape = pvector[1]
  wscale = pvector[2]
  sigmau = pvector[3]
  xi = pvector[4]

  nllh = -lweibullgpd(x, wshape, wscale, u, sigmau, xi, phiu) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}
