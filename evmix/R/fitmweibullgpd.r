#' @export
#' 
#' @title MLE Fitting of Weibull Bulk and GPD Tail Interval Transition Mixture Model
#'
#' @description Maximum likelihood estimation for fitting the extreme valeu 
#' mixture model with the Weibull bulk and GPD tail interval transition mixture model.
#' With options for profile likelihood estimation for threshold and interval half-width,
#' which can both be fixed.
#'
#' @inheritParams fitmnormgpd
#' @inheritParams fweibullgpd
#' @inheritParams fnormgpd
#' @inheritParams fgpd
#' @inheritParams ditmweibullgpd
#' 
#' @details The extreme value mixture model with the Weibull bulk and GPD tail with interval
#' transition is fitted to the entire dataset using maximum likelihood estimation.
#' The estimated parameters, variance-covariance matrix and their standard errors are automatically
#' output.
#' 
#' See \code{\link[evmix:itmweibullgpd]{ditmweibullgpd}} for explanation of Weibull-GPD interval
#' transition model, including mixing functions.
#' 
#' See also help for \code{\link[evmix:fnormgpd]{fnormgpd}} for mixture model fitting details.
#' Only the different features are outlined below for brevity.
#' 
#' The full parameter vector is
#' (\code{wshape}, \code{wscale}, \code{epsilon}, \code{u}, \code{sigmau}, \code{xi})
#' if threshold and interval half-width are both estimated and
#' (\code{wshape}, \code{wscale}, \code{sigmau}, \code{xi})
#' for profile likelihood or fixed threshold and epsilon approach.
#' 
#' If the profile likelihood approach is used, then it is applied to both the threshold and
#' epsilon parameters together. A grid search over all combinations of epsilons and thresholds
#' are considered. The combinations which lead to less than 5 on either side of the interval are 
#' not considered.
#' 
#' A fixed threshold and epsilon approach is acheived by setting a single scalar value to each in 
#' \code{useq} and \code{eseq} respectively.
#' 
#' If the profile likelihood approach is used, then a grid search over all combinations of epsilon and threshold
#' are carried out. The combinations which lead to less than 5 in any any interval are not considered.
#' 
#' Negative data are ignored.
#' 
#' @return Log-likelihood is given by \code{\link[evmix:fitmweibullgpd]{litmweibullgpd}} and it's
#'   wrappers for negative log-likelihood from \code{\link[evmix:fitmweibullgpd]{nlitmweibullgpd}}
#'   and \code{\link[evmix:fitmweibullgpd]{nluitmweibullgpd}}. Profile likelihood for
#'   threshold and interval half-width given by \code{\link[evmix:fitmweibullgpd]{profluitmweibullgpd}}.
#'   Fitting function \code{\link[evmix:fitmweibullgpd]{fitmweibullgpd}} returns a simple list
#'   with the following elements
#'
#' \tabular{ll}{
#'  \code{call}:      \tab \code{optim} call\cr
#'  \code{x}:         \tab data vector \code{x}\cr
#'  \code{init}:      \tab \code{pvector}\cr
#'  \code{fixedeu}:   \tab fixed epsilon and threshold, logical\cr
#'  \code{useq}:      \tab threshold vector for profile likelihood or scalar for fixed threshold\cr
#'  \code{eseq}:      \tab epsilon vector for profile likelihood or scalar for fixed epsilon\cr
#'  \code{nllheuseq}: \tab profile negative log-likelihood at each combination in (eseq, useq)\cr
#'  \code{optim}:     \tab complete \code{optim} output\cr
#'  \code{mle}:       \tab vector of MLE of parameters\cr
#'  \code{cov}:       \tab variance-covariance matrix of MLE of parameters\cr
#'  \code{se}:        \tab vector of standard errors of MLE of parameters\cr
#'  \code{nllh}:      \tab minimum negative log-likelihood\cr
#'  \code{n}:         \tab total sample size\cr
#'  \code{wshape}:    \tab MLE of Weibull shape\cr
#'  \code{wscale}:    \tab MLE of Weibull scale\cr
#'  \code{epsilon}:   \tab MLE of transition half-width\cr
#'  \code{u}:         \tab threshold (fixed or MLE)\cr
#'  \code{sigmau}:    \tab MLE of GPD scale\cr
#'  \code{xi}:        \tab MLE of GPD shape\cr
#' }
#' 
#' @note When \code{pvector=NULL} then the initial values are:
#' \itemize{
#'  \item MLE of Weibull parameters assuming entire population is Weibull; and
#'  \item epsilon is MLE of Weibull standard deviation;
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
#' Holden, L. and Haug, O. (2013). A mixture model for unsupervised tail
#' estimation. arxiv:0902.4137
#' 
#' @author Alfadino Akbar and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @section Acknowledgments: See Acknowledgments in
#'   \code{\link[evmix:fnormgpd]{fnormgpd}}, type \code{help fnormgpd}.
#' 
#' @seealso \code{\link[stats:Weibull]{dweibull}},
#'  \code{\link[evmix:fgpd]{fgpd}} and \code{\link[evmix:gpd]{gpd}}
#'  
#' @aliases fitmweibullgpd litmweibullgpd nlitmweibullgpd profluitmweibullgpd nluitmweibullgpd
#' @family  itmweibullgpd fitmweibullgpd normgpd fnormgpd
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' par(mfrow = c(1, 1))
#' 
#' x = rweibull(1000, shape = 1, scale = 2)
#' xx = seq(-0.2, 10, 0.01)
#' y = dweibull(xx, shape = 1, scale = 2)
#' 
#' # MLE for complete parameter set
#' fit = fitmweibullgpd(x)
#' hist(x, breaks = seq(0, 20, 0.1), freq = FALSE, xlim = c(-0.2, 10))
#' lines(xx, y)
#' with(fit, lines(xx, ditmweibullgpd(xx, wshape, wscale, epsilon, u, sigmau, xi), col="red"))
#' abline(v = fit$u + fit$epsilon * seq(-1, 1), col = "red")
#'   
#' # Profile likelihood for threshold which is then fixed
#' fitfix = fitmweibullgpd(x, eseq = seq(0, 2, 0.1), useq = seq(0.5, 4, 0.1), fixedeu = TRUE)
#' with(fitfix, lines(xx, ditmweibullgpd(xx, wshape, wscale, epsilon, u, sigmau, xi), col="blue"))
#' abline(v = fitfix$u + fitfix$epsilon * seq(-1, 1), col = "blue")
#' legend("topright", c("True Density", "Weibull-GPD ITM", "Profile likelihood"),
#'   col=c("black", "red", "blue"), lty = 1)
#' }
#'   

# maximum likelihood fitting for Weibull bulk with GPD tail
# interval transition mixture model
fitmweibullgpd <- function(x, eseq = NULL, useq = NULL, fixedeu = FALSE, pvector = NULL,
  std.err = TRUE, method = "BFGS", control = list(maxit = 10000), finitelik = TRUE, ...) {

  call <- match.call()
    
  np = 6 # maximum number of parameters

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.posparam(eseq, allowvec = TRUE, allownull = TRUE, allowzero = TRUE)
  check.posparam(useq, allowvec = TRUE, allownull = TRUE)
  check.logic(fixedeu)
  check.logic(std.err)
  check.optim(method)
  check.control(control)
  check.logic(finitelik)

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

  if ((method == "L-BFGS-B") | (method == "BFGS")) finitelik = TRUE
  
  # useq must be specified if epsilon is fixed
  if (fixedeu & is.null(eseq))
    stop("for fixed epsilon approach, eseq must be specified (as scalar or vector)")

  # useq must be specified if threshold is fixed
  if (fixedeu & is.null(useq))
    stop("for fixed threshold approach, useq must be specified (as scalar or vector)")
  
  # Check if profile likelihood or fixed threshold is being used
  # and determine initial values for parameters in each case
  if (is.null(useq)) { # not profile or fixed

    check.nparam(pvector, nparam = np, allownull = TRUE)
    
    if (is.null(pvector)) {
      initfweibull = fitdistr(x, "weibull", lower = c(1e-8, 1e-8))
      pvector[1] = initfweibull$estimate[1]
      pvector[2] = initfweibull$estimate[2]
      pvector[3] = sqrt(pvector[2]^2 * gamma(1 + 2/pvector[1]) - (pvector[2] * gamma(1 + 1/pvector[1]))^2)
      pvector[4] = as.vector(quantile(x, 0.9))
      initfgpd = fgpd(x, pvector[4], std.err = FALSE)
      pvector[5] = initfgpd$sigmau
      pvector[6] = initfgpd$xi
    }
    
  } else { # profile or fixed
    
    check.nparam(pvector, nparam = np - 2, allownull = TRUE)

    # profile likelihood for epsilon and threshold (gridded) or scalar(s) given
    if ((length(useq)*length(eseq)) != 1) {
      
      eugrid = expand.grid(eseq, useq)

      # remove combinations where interval is beyond range of data (must be at least 5 on either side)
      eugrid = eugrid[sapply(eugrid[, 2] + eugrid[, 1], FUN = function(u, x) sum(x > u) >= 5, x = x),]
      eugrid = eugrid[sapply(eugrid[, 2] - eugrid[, 1], FUN = function(u, x) sum(x < u) >= 5, x = x),]
      
      check.posparam(eugrid[, 1], allowvec = TRUE, allowzero = TRUE)
      check.posparam(eugrid[, 2], allowvec = TRUE)
      
      nllheu = apply(eugrid, 1, profleuitmweibullgpd, pvector = pvector, x = x,
        method = method, control = control, finitelik = finitelik, ...)
      
      if (all(!is.finite(nllheu))) stop("thresholds and epsilon combinations are all invalid")

      epsilon = eugrid[which.min(nllheu), 1]
      u = eugrid[which.min(nllheu), 2]

    } else {
      u = useq
      epsilon = eseq
    }

    if (fixedeu) { # threshold and epsilon fixed
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
        pvector[3] = epsilon
        pvector[4] = u
        initfgpd = fgpd(x, pvector[4], std.err = FALSE)
        pvector[5] = initfgpd$sigmau
        pvector[6] = initfgpd$xi
      } else {
        pvector[6] = pvector[4] # shift GPD scale and shape to add in u and epsilon
        pvector[5] = pvector[3]
        pvector[4] = u
        pvector[3] = epsilon
      }
    }
  }

  if (fixedeu) { # fixed threshold and epsilon likelihood
    nllh = nleuitmweibullgpd(pvector, epsilon, u, x)
    
    if (is.infinite(nllh)) {
      pvector[4] = 0.1
      nllh = nleuitmweibullgpd(pvector, epsilon, u, x)
    }
    
    if (is.infinite(nllh)) stop("initial parameter values are invalid")
  
    fit = optim(par = as.vector(pvector), fn = nleuitmweibullgpd, epsilon = epsilon, u = u, x = x,
      finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)    
    
    wshape = fit$par[1]
    wscale = fit$par[2]
    sigmau = fit$par[3]
    xi = fit$par[4]
    
  } else { # complete (non-separable) likelihood
    
    nllh = nlitmweibullgpd(pvector, x)
    
    if (is.infinite(nllh)) {
      pvector[6] = 0.1
      nllh = nlitmweibullgpd(pvector, x)
    }
    
    if (is.infinite(nllh)) stop("initial parameter values are invalid")
  
    fit = optim(par = as.vector(pvector), fn = nlitmweibullgpd, x = x,
      finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)    
    
    wshape = fit$par[1]
    wscale = fit$par[2]
    epsilon = fit$par[3]
    u = fit$par[4]
    sigmau = fit$par[5]
    xi = fit$par[6]
  }
  
  conv = TRUE
  if ((fit$convergence != 0) | any(fit$par == pvector) | (abs(fit$value) >= 1e6)) {
    conv = FALSE
    warning("check convergence")
  }
  
  kappa = 1/(1 + pweibull(u, wshape, wscale))
  
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
    init = as.vector(pvector), fixedeu = fixedeu, useq = useq, eseq = eseq, nllheuseq = nllheu,
    optim = fit, conv = conv, cov = invhess, mle = fit$par, se = se,
    nllh = fit$value, n = n,
    wshape = wshape, wscale = wscale, epsilon = epsilon, u = u, sigmau = sigmau, xi = xi, kappa = kappa)
}

#' @export
#' @aliases fitmweibullgpd litmweibullgpd nlitmweibullgpd profluitmweibullgpd nluitmweibullgpd
#' @rdname  fitmweibullgpd

# log-likelihood function for Weibull bulk with GPD tail
# interval transition mixture model
litmweibullgpd <- function(x, wshape = 1, wscale = 1,
  epsilon = sqrt(wscale^2 * gamma(1 + 2/wshape) - (wscale * gamma(1 + 1/wshape))^2),
  u = qweibull(0.9, wshape, wscale),
  sigmau = sqrt(wscale^2 * gamma(1 + 2/wshape) - (wscale * gamma(1 + 1/wshape))^2),
  xi = 0, log = TRUE) {

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(wshape)
  check.param(wscale)
  check.param(epsilon)
  check.param(u)
  check.param(sigmau)
  check.param(xi)
  check.logic(log)

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

  check.inputn(c(length(wshape), length(wscale), length(epsilon), length(u), length(sigmau), length(xi)),
               allowscalar = TRUE)

  # assume NA or NaN are irrelevant as entire lower tail is now modelled
  # inconsistent with evd library definition
  # hence use which() to ignore these

  xb = x[which(x <= (u - epsilon))]
  xit = x[which((x >= (u - epsilon)) & (x <= (u + epsilon)))]
  xu = x[which(x > (u + epsilon))]
  nb = length(xb)
  nit = length(xit)
  nu = length(xu)

  if ((wscale <= 0) | (wshape <= 0) | (epsilon < 0) | (sigmau <= 0) | (u <= 0) | 
        ((u - epsilon) <= min(x)) | ((u + epsilon) >= max(x))) {
    l = -Inf
  } else {
    kappa = 1/(1 + pweibull(u, wshape, wscale))
    
    if (n != (nb + nit + nu)) {
      stop("total non-finite sample size is not equal to those above/below interval or within it")    
    }

    syu = 1 + xi * (xu - u) / sigmau
  
    if (min(syu) <= 0) {
      l = -Inf
    } else { 
       l = lgpd(xu, u, sigmau, xi)
       l = l + (wshape - 1) * sum(log(xb)) - sum(xb^wshape) / wscale^wshape + nb * log(wshape) - nb * wshape * log(wscale)
       if (nit > 0) l = l + sum(ditmweibullgpd(xit, wshape, wscale, epsilon, u, sigmau, xi, log = TRUE))
       l = l + (n - nit) * log(kappa)
    }
  }
  
  if (!log) l = exp(l)
  
  l
}

#' @export
#' @aliases fitmweibullgpd litmweibullgpd nlitmweibullgpd profluitmweibullgpd nluitmweibullgpd
#' @rdname  fitmweibullgpd

# negative log-likelihood function for Weibull bulk with GPD tail
# interval transition mixture model
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlitmweibullgpd <- function(pvector, x, finitelik = FALSE) {

  np = 6 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.logic(finitelik)

  wshape = pvector[1]
  wscale = pvector[2]
  epsilon = pvector[3]
  u = pvector[4]
  sigmau = pvector[5]
  xi = pvector[6]

  nllh = -litmweibullgpd(x, wshape, wscale, epsilon, u, sigmau, xi) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}

#' @export
#' @aliases fitmweibullgpd litmweibullgpd nlitmweibullgpd profluitmweibullgpd nluitmweibullgpd
#' @rdname  fitmweibullgpd

# profile negative log-likelihood function for given threshold and epsilon for
# Weibull bulk with GPD tail interval transition mixture model
# designed for sapply to loop over matrix with two columns vector of threshold and epsilon pairs
# (hence eu is first input)
profleuitmweibullgpd <- function(eu, pvector, x, method = "BFGS",
  control = list(maxit = 10000), finitelik = TRUE, ...) {

  np = 6 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np - 2, allownull = TRUE)
  check.nparam(eu, nparam = 2);
  check.posparam(eu[1], allowzero = TRUE)
  check.posparam(eu[2])
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.optim(method)
  check.control(control)
  check.logic(finitelik)

  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    x = x[is.finite(x)] # ignore missing and infinite cases
  }

  if (any(x < 0)) {
    warning("negative values have been removed")
    x = x[x >= 0]
  }

  check.quant(x)

  # check initial values for other parameters, try usual alternative
  if (!is.null(pvector)) {
    nllh = nleuitmweibullgpd(pvector, eu[1], eu[2], x)
    
    if (is.infinite(nllh)) pvector = NULL
  }

  if (is.null(pvector)) {
    initfweibull = fitdistr(x, "weibull", lower = c(1e-8, 1e-8))
    pvector[1] = initfweibull$estimate[1]
    pvector[2] = initfweibull$estimate[2]
    initfgpd = fgpd(x, eu[2], std.err = FALSE)
    pvector[3] = initfgpd$sigmau
    pvector[4] = initfgpd$xi
    nllh = nleuitmweibullgpd(pvector, eu[1], eu[2], x)
  }
    
  if (is.infinite(nllh)) {
    pvector[4] = 0.1
    nllh = nleuitmweibullgpd(pvector, eu[1], eu[2], x)
  }
  
  # if still invalid then output cleanly
  if (is.infinite(nllh)) {
    warning(paste("initial parameter values for threshold u =", eu[2], "and epsilon=", eu[1], "are invalid"))
    fit = list(par = rep(NA, np), value = Inf, counts = 0, convergence = NA, 
      message = "initial values invalid", hessian = rep(NA, np))
  } else {

    fit = optim(par = as.vector(pvector), fn = nleuitmweibullgpd, epsilon = eu[1], u = eu[2], x = x,
      finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)
  }
    
  if (finitelik & is.infinite(fit$value)) {
    fit$value = sign(fit$value) * 1e6
  }

  fit$value
}

#' @export
#' @aliases fitmweibullgpd litmweibullgpd nlitmweibullgpd profluitmweibullgpd nluitmweibullgpd
#' @rdname  fitmweibullgpd

# negative log-likelihood function for Weibull bulk with GPD tail
# interval transition mixture model
# (wrapper for likelihood, designed for threshold to be fixed and other parameters optimised)
nleuitmweibullgpd <- function(pvector, epsilon, u, x, finitelik = FALSE) {

  np = 6 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np - 2)
  check.posparam(epsilon, allowzero = TRUE)
  check.posparam(u)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.logic(finitelik)
    
  wshape = pvector[1]
  wscale = pvector[2]
  sigmau = pvector[3]
  xi = pvector[4]

  nllh = -litmweibullgpd(x, wshape, wscale, epsilon, u, sigmau, xi)
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}
