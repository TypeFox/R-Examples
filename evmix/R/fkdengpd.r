#' @export
#' 
#' @title MLE Fitting of Kernel Density Estimate for Bulk and GPD Tail Extreme Value Mixture Model
#'
#' @description Maximum likelihood estimation for fitting the extreme value 
#' mixture model with kernel density estimate for bulk distribution upto the threshold and conditional
#' GPD above threshold. With options for profile likelihood estimation for threshold and
#' fixed threshold approach.
#'
#' @param bw           scalar bandwidth for kernel (as standard deviations of kernel)
#' @param lambda       scalar bandwidth for kernel (as half-width of kernel)
#' @inheritParams fnormgpd
#' @inheritParams fgpd
#' @inheritParams fkden
#' 
#' @details The extreme value mixture model with kernel density estimate for bulk and GPD tail is 
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
#' Cross-validation likelihood is used for KDE, but standard likelihood is used
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
#' @section Warning:
#' See important warnings about cross-validation likelihood estimation in 
#' \code{\link[evmix:fkden]{fkden}}, type \code{help fkden}.
#' 
#' @return Log-likelihood is given by \code{\link[evmix:fkdengpd]{lkdengpd}} and it's
#'   wrappers for negative log-likelihood from \code{\link[evmix:fkdengpd]{nlkdengpd}}
#'   and \code{\link[evmix:fkdengpd]{nlukdengpd}}. Profile likelihood for single
#'   threshold given by \code{\link[evmix:fkdengpd]{proflukdengpd}}. Fitting function
#'   \code{\link[evmix:fkdengpd]{fkdengpd}} returns a simple list with the
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
#'  \code{lambda}:    \tab MLE of lambda (kernel half-width)\cr
#'  \code{u}:         \tab threshold (fixed or MLE)\cr
#'  \code{sigmau}:    \tab MLE of GPD scale\cr
#'  \code{xi}:        \tab MLE of GPD shape\cr
#'  \code{phiu}:      \tab MLE of tail fraction (bulk model or parameterised approach)\cr
#'  \code{se.phiu}:   \tab standard error of MLE of tail fraction\cr
#'  \code{bw}:        \tab MLE of bw (kernel standard deviations)\cr
#'  \code{kernel}:    \tab kernel name\cr
#' }
#' 
#' @note The data and kernel centres are both vectors. Infinite and missing sample values
#' (and kernel centres) are dropped.
#' 
#' When \code{pvector=NULL} then the initial values are:
#' \itemize{
#'  \item normal reference rule for bandwidth, using the \code{\link[stats:bandwidth]{bw.nrd0}} function, which is
#'    consistent with the \code{\link[stats:density]{density}} function. At least two kernel
#'    centres must be provided as the variance needs to be estimated.
#'  \item threshold 90\% quantile (not relevant for profile likelihood for threshold or fixed threshold approaches);
#'  \item MLE of GPD parameters above threshold. 
#' }
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
#'  
#' @aliases fkdengpd lkdengpd nlkdengpd proflukdengpd nlukdengpd
#' @family  kdengpd kdengpdcon fkdengpd fkdengpdcon normgpd fnormgpd
#'          kden bckden bckdengpd bckdengpdcon fkden fbckden fbckdengpd fbckdengpdcon
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' par(mfrow = c(2, 1))
#' 
#' x = rnorm(1000)
#' xx = seq(-4, 4, 0.01)
#' y = dnorm(xx)
#' 
#' # Bulk model based tail fraction
#' fit = fkdengpd(x)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-4, 4))
#' lines(xx, y)
#' with(fit, lines(xx, dkdengpd(xx, x, lambda, u, sigmau, xi), col="red"))
#' abline(v = fit$u, col = "red")
#'   
#' # Parameterised tail fraction
#' fit2 = fkdengpd(x, phiu = FALSE)
#' with(fit2, lines(xx, dkdengpd(xx, x, lambda, u, sigmau, xi, phiu), col="blue"))
#' abline(v = fit2$u, col = "blue")
#' legend("topright", c("True Density","Bulk Tail Fraction","Parameterised Tail Fraction"),
#'   col=c("black", "red", "blue"), lty = 1)
#'   
#' # Profile likelihood for initial value of threshold and fixed threshold approach
#' fitu = fkdengpd(x, useq = seq(0, 2, length = 20))
#' fitfix = fkdengpd(x, useq = seq(0, 2, length = 20), fixedu = TRUE)
#' 
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-4, 4))
#' lines(xx, y)
#' with(fit, lines(xx, dkdengpd(xx, x, lambda, u, sigmau, xi), col="red"))
#' abline(v = fit$u, col = "red")
#' with(fitu, lines(xx, dkdengpd(xx, x, lambda, u, sigmau, xi), col="purple"))
#' abline(v = fitu$u, col = "purple")
#' with(fitfix, lines(xx, dkdengpd(xx, x, lambda, u, sigmau, xi), col="darkgreen"))
#' abline(v = fitfix$u, col = "darkgreen")
#' legend("topright", c("True Density","Default initial value (90% quantile)",
#'  "Prof. lik. for initial value", "Prof. lik. for fixed threshold"),
#'  col=c("black", "red", "purple", "darkgreen"), lty = 1)
#' }
#'   

# maximum likelihood fitting for kernel density estimate for bulk with GPD for upper tail
fkdengpd <- function(x, phiu = TRUE, useq = NULL, fixedu = FALSE, pvector = NULL, kernel = "gaussian",
  add.jitter = FALSE, factor = 0.1, amount = NULL,
  std.err = TRUE, method = "BFGS", control = list(maxit = 10000), finitelik = TRUE, ...) {

  call <- match.call()
    
  np = 4 # maximum number of parameters

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.logic(phiu)
  check.param(useq, allowvec = TRUE, allownull = TRUE)
  check.logic(fixedu)
  check.logic(std.err)
  check.optim(method)
  check.control(control)
  check.logic(finitelik)

  check.kernel(kernel)
  check.posparam(factor)
  check.posparam(amount, allownull = TRUE)
  check.logic(add.jitter)
  
  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    x = x[is.finite(x)] # ignore missing and infinite cases
  }

  check.quant(x)
  n = length(x)

  if (add.jitter) x = jitter(x, factor, amount)

  xuniq = unique(x)
  if (length(xuniq) < (0.95*n))
    warning("data may be rounded, as more than 5% are ties, so bandwidth could be biased to zero")

  if ((method == "L-BFGS-B") | (method == "BFGS")) finitelik = TRUE
  
  # useq must be specified if threshold is fixed
  if (fixedu & is.null(useq))
    stop("for fixed threshold approach, useq must be specified (as scalar or vector)")
  
  # Check if profile likelihood or fixed threshold is being used
  # and determine initial values for parameters in each case
  if (is.null(useq)) { # not profile or fixed

    check.nparam(pvector, nparam = np, allownull = TRUE)
    
    if (is.null(pvector)) {
      if (n == 1) {
        stop("Automated bandwidth estimation requires 2 or more kernel centres")
      } else if (n < 10) {
        stop("Automated bandwidth estimation unreliable with less than 10 kernel centres")
      }
      pvector[1] = klambda(bw.nrd0(x), kernel)
      pvector[2] = as.vector(quantile(x, 0.9))
      initfgpd = fgpd(x, pvector[2], std.err = FALSE)
      pvector[3] = initfgpd$sigmau
      pvector[4] = initfgpd$xi
    }
    
  } else { # profile or fixed
    
    check.nparam(pvector, nparam = np - 1, allownull = TRUE)

    # profile likelihood for threshold or scalar given
    if (length(useq) != 1) {
      
      # remove thresholds with less than 5 excesses
      useq = useq[sapply(useq, FUN = function(u, x) sum(x > u) > 5, x = x)]
      check.param(useq, allowvec = TRUE)
      
      nllhu = sapply(useq, proflukdengpd, pvector = pvector, x = x, phiu = phiu, kernel = kernel,
        method = method, control = control, finitelik = finitelik, ...)
      
      if (all(!is.finite(nllhu))) stop("thresholds are all invalid")
      u = useq[which.min(nllhu)]

    } else {
      u = useq
    }

    if (fixedu) { # threshold fixed
      if (is.null(pvector)) {
        if (n == 1) {
          stop("Automated bandwidth estimation requires 2 or more kernel centres")
        } else if (n < 10) {
          stop("Automated bandwidth estimation unreliable with less than 10 kernel centres")
        }
        pvector[1] = klambda(bw.nrd0(x), kernel)
        initfgpd = fgpd(x, u, std.err = FALSE)
        pvector[2] = initfgpd$sigmau
        pvector[3] = initfgpd$xi
      }
    } else { # threshold as initial value in usual MLE
      if (is.null(pvector)) {
        if (n == 1) {
          stop("Automated bandwidth estimation requires 2 or more kernel centres")
        } else if (n < 10) {
          stop("Automated bandwidth estimation unreliable with less than 10 kernel centres")
        }
        pvector[1] = klambda(bw.nrd0(x), kernel)
        pvector[2] = u
        initfgpd = fgpd(x, pvector[2], std.err = FALSE)
        pvector[3] = initfgpd$sigmau
        pvector[4] = initfgpd$xi
      } else {
        pvector[4] = pvector[3] # shift GPD scale and shape to add in u
        pvector[3] = pvector[2]
        pvector[2] = u
      }
    }
  }

  if (fixedu) { # fixed threshold (separable) likelihood
    nllh = nlukdengpd(pvector, u, x, phiu, kernel = kernel)
    if (is.infinite(nllh)) {
      pvector[3] = 0.1
      nllh = nlukdengpd(pvector, u, x, phiu, kernel = kernel)    
    }
    
    if (is.infinite(nllh)) stop("initial parameter values are invalid")
  
    fit = optim(par = as.vector(pvector), fn = nlukdengpd, u = u, x = x, phiu = phiu, kernel = kernel,
      finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)    
    
    lambda = fit$par[1]
    sigmau = fit$par[2]
    xi = fit$par[3]
    
  } else { # complete (non-separable) likelihood
    
    nllh = nlkdengpd(pvector, x, phiu, kernel = kernel)
    if (is.infinite(nllh)) {
      pvector[4] = 0.1
      nllh = nlkdengpd(pvector, x, phiu, kernel = kernel)    
    }
    if (is.infinite(nllh)) stop("initial parameter values are invalid")
  
    fit = optim(par = as.vector(pvector), fn = nlkdengpd, x = x, phiu = phiu, kernel = kernel,
      finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)    
    
    lambda = fit$par[1]
    u = fit$par[2]
    sigmau = fit$par[3]
    xi = fit$par[4]
  }

  bw = kbw(fit$par[1], kernel)

  conv = TRUE
  if ((fit$convergence != 0) | any(fit$par == pvector) | (abs(fit$value) >= 1e6)) {
    conv = FALSE
    warning("check convergence")
  }

  pu = pkdenx(u, x, lambda, kernel)
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
    lambda = lambda, u = u, sigmau = sigmau, xi = xi, phiu = phiu, se.phiu = se.phiu, bw = bw, kernel = kernel)
}

#' @export
#' @aliases fkdengpd lkdengpd nlkdengpd proflukdengpd nlukdengpd
#' @rdname  fkdengpd

# log-likelihood function for kernel density estimate for bulk with GPD for upper tail
# cross-validation for KDE component
lkdengpd <- function(x, lambda = NULL, u = 0, sigmau = 1, xi = 0, phiu = TRUE,
  bw = NULL, kernel = "gaussian", log = TRUE) {

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

  kernel = ifelse(kernel == "rectangular", "uniform", kernel)
  kernel = ifelse(kernel == "normal", "gaussian", kernel)
  
  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    x = x[is.finite(x)] # ignore missing and infinite cases
  }

  check.quant(x)
  n = length(x)

  if (is.null(lambda) & is.null(bw)) {
    if (n == 1) {
      stop("Automated bandwidth estimation requires 2 or more kernel centres")
    } else if (n < 10) {
      warning("Automated bandwidth estimation unreliable with less than 10 kernel centres")
    }
    bw = bw.nrd0(x)
  }

  if (is.null(lambda)) lambda = klambda(bw, kernel, lambda)
  
  check.inputn(c(length(lambda), length(u), length(sigmau), length(xi), length(phiu)), allowscalar = TRUE)

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

  if ((lambda <= 0) | (sigmau <= 0) | (u <= min(x)) | (u >= max(x))) {
    l = -Inf
  } else {
    if (is.logical(phiu)) {
      pu = pkdenx(u, x, lambda, kernel)
      if (phiu) {
        phiu = 1 - pu
      } else {
        phiu = nu / n
      }
    }
    phib = (1 - phiu) / pu
  
    syu = 1 + xi * (xu - u) / sigmau  
  
    if ((min(syu) <= 0) | (phiu <= 0) | (phiu >= 1) | (pu <= 0) | (pu >= 1)) {
      l = -Inf
    } else { 
      l = lgpd(xu, u, sigmau, xi, phiu)
      l = l + lkden(xb, lambda, kernel = kernel, extracentres = xu) + nb*log(phib)
    }
  }
  
  if (!log) l = exp(l)
  
  l
}

#' @export
#' @aliases fkdengpd lkdengpd nlkdengpd proflukdengpd nlukdengpd
#' @rdname  fkdengpd

# negative log-likelihood function for kernel density estimate for bulk with GPD for upper tail
# cross-validation for KDE component
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlkdengpd <- function(pvector, x, phiu = TRUE, kernel = "gaussian", finitelik = FALSE) {

  np = 4 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.phiu(phiu, allowfalse = TRUE)
  check.logic(finitelik)

  check.kernel(kernel)
  kernel = ifelse(kernel == "rectangular", "uniform", kernel)
  kernel = ifelse(kernel == "normal", "gaussian", kernel)

  lambda = pvector[1]
  u = pvector[2]
  sigmau = pvector[3]
  xi = pvector[4]

  nllh = -lkdengpd(x, lambda, u, sigmau, xi, phiu, kernel = kernel) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}

#' @export
#' @aliases fkdengpd lkdengpd nlkdengpd proflukdengpd nlukdengpd
#' @rdname  fkdengpd

# profile negative log-likelihood function for given threshold for
# kernel density estimate for bulk with GPD for upper tail
# designed for sapply to loop over vector of thresholds (hence u is first input)
# cross-validation for KDE component
proflukdengpd <- function(u, pvector, x, phiu = TRUE, kernel = "gaussian",
  method = "BFGS", control = list(maxit = 10000), finitelik = TRUE, ...) {

  np = 4 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np - 1, allownull = TRUE)
  check.param(u)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.phiu(phiu, allowfalse = TRUE)
  check.optim(method)
  check.control(control)
  check.logic(finitelik)

  check.kernel(kernel)

  kernel = ifelse(kernel == "rectangular", "uniform", kernel)
  kernel = ifelse(kernel == "normal", "gaussian", kernel)

  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    x = x[is.finite(x)] # ignore missing and infinite cases
  }

  check.quant(x)
  n = length(x)
  
  # check initial values for other parameters, try usual alternative
  if (!is.null(pvector)) {
    nllh = nlukdengpd(pvector, u, x, phiu, kernel = kernel)
    
    if (is.infinite(nllh)) pvector = NULL
  }

  if (is.null(pvector)) {
    if (n == 1) {
      stop("Automated bandwidth estimation requires 2 or more kernel centres")
    } else if (n < 10) {
        stop("Automated bandwidth estimation unreliable with less than 10 kernel centres")
    }
    pvector[1] = klambda(bw.nrd0(x), kernel)
    initfgpd = fgpd(x, u, std.err = FALSE)
    pvector[2] = initfgpd$sigmau
    pvector[3] = initfgpd$xi
    nllh = nlukdengpd(pvector, u, x, phiu, kernel = kernel)
  }

  if (is.infinite(nllh)) {
    pvector[3] = 0.1
    nllh = nlukdengpd(pvector, u, x, phiu, kernel = kernel)    
  }
  
  # if still invalid then output cleanly
  if (is.infinite(nllh)) {
    warning(paste("initial parameter values for threshold u =", u, "are invalid"))
    fit = list(par = rep(NA, np), value = Inf, counts = 0, convergence = NA, 
      message = "initial values invalid", hessian = rep(NA, np))
  } else {

    fit = optim(par = as.vector(pvector), fn = nlukdengpd, u = u, x = x, phiu = phiu, kernel = kernel,
    finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)
  }
    
  if (finitelik & is.infinite(fit$value)) {
    fit$value = sign(fit$value) * 1e6
  }

  fit$value
}

#' @export
#' @aliases fkdengpd lkdengpd nlkdengpd proflukdengpd nlukdengpd
#' @rdname  fkdengpd

# negative log-likelihood function for kernel density estimate for bulk with GPD for upper tail
# (wrapper for likelihood, designed for threshold to be fixed and other parameters optimised)
# cross-validation for KDE component
nlukdengpd <- function(pvector, u, x, phiu = TRUE, kernel = "gaussian", finitelik = FALSE) {

  np = 4 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np - 1)
  check.param(u)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.phiu(phiu, allowfalse = TRUE)
  check.logic(finitelik)

  check.kernel(kernel)
  kernel = ifelse(kernel == "rectangular", "uniform", kernel)
  kernel = ifelse(kernel == "normal", "gaussian", kernel)
    
  lambda = pvector[1]
  sigmau = pvector[2]
  xi = pvector[3]

  nllh = -lkdengpd(x, lambda, u, sigmau, xi, phiu, kernel = kernel) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}
