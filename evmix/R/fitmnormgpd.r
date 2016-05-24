#' @export
#' 
#' @title MLE Fitting of Normal Bulk and GPD Tail Interval Transition Mixture Model
#'
#' @description Maximum likelihood estimation for fitting the extreme value 
#' mixture model with the normal bulk and GPD tail interval transition mixture model.
#' With options for profile likelihood estimation for threshold and interval half-width,
#' which can both be fixed.
#'
#' @param useq    vector of thresholds (or scalar) to be considered in profile likelihood or
#'                \code{NULL} for no profile likelihood
#' @param eseq    vector of epsilons (or scalar) to be considered in profile likelihood or
#'                \code{NULL} for no profile likelihood
#' @param fixedeu logical, should threshold and epsilon be fixed
#'                (at either scalar value in \code{useq} and \code{eseq},
#'                or estimated from maximum of profile likelihood evaluated at
#'                grid of thresholds and epsilons in \code{useq} and \code{eseq})
#' @param eu      vector of epsilon and threshold pair considered in profile likelihood
#' @inheritParams fnormgpd
#' @inheritParams fnormgpd
#' @inheritParams fgpd
#' @inheritParams ditmnormgpd
#' 
#' @details The extreme value mixture model with the normal bulk and GPD tail with interval
#' transition is fitted to the entire dataset using maximum likelihood estimation.
#' The estimated parameters, variance-covariance matrix and their standard errors are automatically
#' output.
#' 
#' See \code{\link[evmix:itmnormgpd]{ditmnormgpd}} for explanation of normal-GPD interval
#' transition model, including mixing functions.
#' 
#' See also help for \code{\link[evmix:fnormgpd]{fnormgpd}} for mixture model fitting details.
#' Only the different features are outlined below for brevity.
#' 
#' The full parameter vector is
#' (\code{nmean}, \code{nsd}, \code{epsilon}, \code{u}, \code{sigmau}, \code{xi})
#' if threshold and interval half-width are both estimated and
#' (\code{nmean}, \code{nsd}, \code{sigmau}, \code{xi})
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
#' @return Log-likelihood is given by \code{\link[evmix:fitmnormgpd]{litmnormgpd}} and it's
#'   wrappers for negative log-likelihood from \code{\link[evmix:fitmnormgpd]{nlitmnormgpd}}
#'   and \code{\link[evmix:fitmnormgpd]{nluitmnormgpd}}. Profile likelihood for
#'   threshold and interval half-width given by \code{\link[evmix:fitmnormgpd]{profluitmnormgpd}}.
#'   Fitting function \code{\link[evmix:fitmnormgpd]{fitmnormgpd}} returns a simple list
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
#'  \code{nmean}:     \tab MLE of normal shape\cr
#'  \code{nsd}:       \tab MLE of normal scale\cr
#'  \code{epsilon}:   \tab MLE of transition half-width\cr
#'  \code{u}:         \tab threshold (fixed or MLE)\cr
#'  \code{sigmau}:    \tab MLE of GPD scale\cr
#'  \code{xi}:        \tab MLE of GPD shape\cr
#' }
#' 
#' @note When \code{pvector=NULL} then the initial values are:
#' \itemize{
#'  \item MLE of normal parameters assuming entire population is normal; and
#'  \item epsilon is MLE of normal standard deviation;
#'  \item threshold 90\% quantile (not relevant for profile likelihood for threshold or fixed threshold approaches);
#'  \item MLE of GPD parameters above threshold. 
#' }
#' 
#' @references
#' \url{http://www.math.canterbury.ac.nz/~c.scarrott/evmix}
#' 
#' \url{http://en.wikipedia.org/wiki/normal_distribution}
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
#' @seealso \code{\link[evmix:fnormgpd]{fnormgpd}}, \code{\link[stats:Normal]{dnorm}},
#'  \code{\link[evmix:fgpd]{fgpd}} and \code{\link[evmix:gpd]{gpd}}
#'  
#' @aliases fitmnormgpd litmnormgpd nlitmnormgpd profluitmnormgpd nluitmnormgpd
#' @family  itmnormgpd fitmnormgpd normgpd fnormgpd
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
#' # MLE for complete parameter set
#' fit = fitmnormgpd(x)
#' hist(x, breaks = seq(-6, 6, 0.1), freq = FALSE, xlim = c(-4, 4))
#' lines(xx, y)
#' with(fit, lines(xx, ditmnormgpd(xx, nmean, nsd, epsilon, u, sigmau, xi), col="red"))
#' abline(v = fit$u + fit$epsilon * seq(-1, 1), col = "red")
#'   
#' # Profile likelihood for threshold which is then fixed
#' fitfix = fitmnormgpd(x, eseq = seq(0, 2, 0.1), useq = seq(0, 2.5, 0.1), fixedeu = TRUE)
#' with(fitfix, lines(xx, ditmnormgpd(xx, nmean, nsd, epsilon, u, sigmau, xi), col="blue"))
#' abline(v = fitfix$u + fitfix$epsilon * seq(-1, 1), col = "blue")
#' legend("topright", c("True Density", "normal-GPD ITM", "Profile likelihood"),
#'   col=c("black", "red", "blue"), lty = 1)
#' }
#'   

# maximum likelihood fitting for normal bulk with GPD tail
# interval transition mixture model
fitmnormgpd <- function(x, eseq = NULL, useq = NULL, fixedeu = FALSE, pvector = NULL,
  std.err = TRUE, method = "BFGS", control = list(maxit = 10000), finitelik = TRUE, ...) {

  call <- match.call()
    
  np = 6 # maximum number of parameters

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.posparam(eseq, allowvec = TRUE, allownull = TRUE, allowzero = TRUE)
  check.param(useq, allowvec = TRUE, allownull = TRUE)
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
  if (fixedeu & is.null(useq))
    stop("for fixed threshold approach, useq must be specified (as scalar or vector)")
  
  # Check if profile likelihood or fixed threshold is being used
  # and determine initial values for parameters in each case
  if (is.null(useq)) { # not profile or fixed

    check.nparam(pvector, nparam = np, allownull = TRUE)
    
    if (is.null(pvector)) {
      pvector[1] = mean(x)
      pvector[2] = sd(x)
      pvector[3] = pvector[2]
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
      check.param(eugrid[, 2], allowvec = TRUE)

      nllheu = apply(eugrid, 1, profleuitmnormgpd, pvector = pvector, x = x,
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
        pvector[1] = mean(x)
        pvector[2] = sd(x)
        initfgpd = fgpd(x, u, std.err = FALSE)
        pvector[3] = initfgpd$sigmau
        pvector[4] = initfgpd$xi
      }
    } else { # threshold as initial value in usual MLE
      if (is.null(pvector)) {
        pvector[1] = mean(x)
        pvector[2] = sd(x)
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
    nllh = nleuitmnormgpd(pvector, epsilon, u, x)
    
    if (is.infinite(nllh)) {
      pvector[4] = 0.1
      nllh = nleuitmnormgpd(pvector, epsilon, u, x)
    }
    
    if (is.infinite(nllh)) stop("initial parameter values are invalid")
  
    fit = optim(par = as.vector(pvector), fn = nleuitmnormgpd, epsilon = epsilon, u = u, x = x,
      finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)    
    
    nmean = fit$par[1]
    nsd = fit$par[2]
    sigmau = fit$par[3]
    xi = fit$par[4]
    
  } else { # complete (non-separable) likelihood
    
    nllh = nlitmnormgpd(pvector, x)
    
    if (is.infinite(nllh)) {
      pvector[6] = 0.1
      nllh = nlitmnormgpd(pvector, x)
    }
    
    if (is.infinite(nllh)) stop("initial parameter values are invalid")
  
    fit = optim(par = as.vector(pvector), fn = nlitmnormgpd, x = x,
      finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)    
    
    nmean = fit$par[1]
    nsd = fit$par[2]
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
  
  kappa = 1/(1 + pnorm(u, nmean, nsd))
  
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
    nmean = nmean, nsd = nsd, epsilon = epsilon, u = u, sigmau = sigmau, xi = xi, kappa = kappa)
}

#' @export
#' @aliases fitmnormgpd litmnormgpd nlitmnormgpd profluitmnormgpd nluitmnormgpd
#' @rdname  fitmnormgpd

# log-likelihood function for normal bulk with GPD tail
# interval transition mixture model
litmnormgpd <- function(x, nmean = 0, nsd = 1, epsilon = nsd,  u = qnorm(0.9, nmean, nsd),
  sigmau = nsd, xi = 0, log = TRUE) {

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(nmean)
  check.param(nsd)
  check.param(epsilon)
  check.param(u)
  check.param(sigmau)
  check.param(xi)
  check.logic(log)

  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    x = x[is.finite(x)] # ignore missing and infinite cases
  }

  check.quant(x)
  n = length(x)

  check.inputn(c(length(nmean), length(nsd), length(epsilon), length(u), length(sigmau), length(xi)),
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

  if ((nsd <= 0) | (epsilon < 0) | (sigmau <= 0) | ((u - epsilon) <= min(x)) | ((u + epsilon) >= max(x))) {
    l = -Inf
  } else {
    kappa = 1/(1 + pnorm(u, nmean, nsd))
    
    if (n != (nb + nit + nu)) {
      stop("total non-finite sample size is not equal to those above/below interval or within it")    
    }

    syu = 1 + xi * (xu - u) / sigmau
    yb = (xb - nmean) / nsd # used for normal
  
    if (min(syu) <= 0) {
      l = -Inf
    } else { 
       l = lgpd(xu, u, sigmau, xi)
       l = l - nb * log(2 * pi * nsd ^ 2) / 2 - sum(yb ^ 2) / 2
       if (nit > 0) l = l + sum(ditmnormgpd(xit, nmean, nsd, epsilon, u, sigmau, xi, log = TRUE))
       l = l + (n - nit) * log(kappa)
    }
  }
  
  if (!log) l = exp(l)
  
  l
}

#' @export
#' @aliases fitmnormgpd litmnormgpd nlitmnormgpd profluitmnormgpd nluitmnormgpd
#' @rdname  fitmnormgpd

# negative log-likelihood function for normal bulk with GPD tail
# interval transition mixture model
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlitmnormgpd <- function(pvector, x, finitelik = FALSE) {

  np = 6 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.logic(finitelik)

  nmean = pvector[1]
  nsd = pvector[2]
  epsilon = pvector[3]
  u = pvector[4]
  sigmau = pvector[5]
  xi = pvector[6]

  nllh = -litmnormgpd(x, nmean, nsd, epsilon, u, sigmau, xi) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}

#' @export
#' @aliases fitmnormgpd litmnormgpd nlitmnormgpd profluitmnormgpd nluitmnormgpd
#' @rdname  fitmnormgpd

# profile negative log-likelihood function for given threshold and epsilon for
# normal bulk with GPD tail interval transition mixture model
# designed for sapply to loop over matrix with two columns vector of threshold and epsilon pairs
# (hence eu is first input)
profleuitmnormgpd <- function(eu, pvector, x, method = "BFGS",
  control = list(maxit = 10000), finitelik = TRUE, ...) {

  np = 6 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np - 2, allownull = TRUE)
  check.nparam(eu, nparam = 2);
  check.posparam(eu[1], allowzero = TRUE)
  check.param(eu[2])
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
    nllh = nleuitmnormgpd(pvector, eu[1], eu[2], x)
    
    if (is.infinite(nllh)) pvector = NULL
  }

  if (is.null(pvector)) {
    pvector[1] = mean(x)
    pvector[2] = sd(x)
    initfgpd = fgpd(x, eu[2], std.err = FALSE)
    pvector[3] = initfgpd$sigmau
    pvector[4] = initfgpd$xi
    nllh = nleuitmnormgpd(pvector, eu[1], eu[2], x)
  }

  if (is.infinite(nllh)) {
    pvector[4] = 0.1
    nllh = nleuitmnormgpd(pvector, eu[1], eu[2], x)
  }

  # if still invalid then output cleanly
  if (is.infinite(nllh)) {
    warning(paste("initial parameter values for threshold u =", eu[2], "and epsilon=", eu[1], "are invalid"))
    fit = list(par = rep(NA, np), value = Inf, counts = 0, convergence = NA, 
      message = "initial values invalid", hessian = rep(NA, np))
  } else {

    fit = optim(par = as.vector(pvector), fn = nleuitmnormgpd, epsilon = eu[1], u = eu[2], x = x,
      finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)
  }
    
  if (finitelik & is.infinite(fit$value)) {
    fit$value = sign(fit$value) * 1e6
  }

  fit$value
}

#' @export
#' @aliases fitmnormgpd litmnormgpd nlitmnormgpd profluitmnormgpd nluitmnormgpd
#' @rdname  fitmnormgpd

# negative log-likelihood function for normal bulk with GPD tail
# interval transition mixture model
# (wrapper for likelihood, designed for threshold to be fixed and other parameters optimised)
nleuitmnormgpd <- function(pvector, epsilon, u, x, finitelik = FALSE) {

  np = 6 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np - 2)
  check.posparam(epsilon, allowzero = TRUE)
  check.param(u)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.logic(finitelik)
    
  nmean = pvector[1]
  nsd = pvector[2]
  sigmau = pvector[3]
  xi = pvector[4]

  nllh = -litmnormgpd(x, nmean, nsd, epsilon, u, sigmau, xi)
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}
