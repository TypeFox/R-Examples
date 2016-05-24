#' @export
#' 
#' @title MLE Fitting of Normal Bulk and GPD Tail Extreme Value Mixture Model with Single Continuity Constraint
#'
#' @description Maximum likelihood estimation for fitting the extreme value 
#' mixture model with normal for bulk distribution upto the threshold and conditional
#' GPD above threshold with continuity at threshold. With options for profile likelihood estimation for threshold and
#' fixed threshold approach.
#'
#' @inheritParams fnormgpd
#' 
#' @details The extreme value mixture model with normal bulk and GPD tail with continuity at threshold is 
#' fitted to the entire dataset using maximum likelihood estimation. The estimated
#' parameters, variance-covariance matrix and their standard errors are automatically
#' output.
#' 
#' See help for \code{\link[evmix:fnormgpd]{fnormgpd}} for full details, type \code{help fnormgpd}. Only
#' the different features are outlined below for brevity.
#' 
#' The GPD \code{sigmau} parameter is now specified as function of other parameters, see 
#' help for \code{\link[evmix:normgpdcon]{dnormgpdcon}} for details, type \code{help normgpdcon}.
#' Therefore, \code{sigmau} should not be included in the parameter vector if initial values
#' are provided, making the full parameter vector 
#' (\code{nmean}, \code{nsd}, \code{u}, \code{xi}) if threshold is also estimated and
#' (\code{nmean}, \code{nsd}, \code{xi}) for profile likelihood or fixed threshold approach.
#' 
#' @return Log-likelihood is given by \code{\link[evmix:fnormgpdcon]{lnormgpdcon}} and it's
#'   wrappers for negative log-likelihood from \code{\link[evmix:fnormgpdcon]{nlnormgpdcon}}
#'   and \code{\link[evmix:fnormgpdcon]{nlunormgpdcon}}. Profile likelihood for single
#'   threshold given by \code{\link[evmix:fnormgpdcon]{proflunormgpdcon}}. Fitting function
#'   \code{\link[evmix:fnormgpdcon]{fnormgpdcon}} returns a simple list with the
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
#'  \code{nmean}:     \tab MLE of normal mean\cr
#'  \code{nsd}:       \tab MLE of normal standard deviation\cr
#'  \code{u}:         \tab threshold (fixed or MLE)\cr
#'  \code{sigmau}:    \tab MLE of GPD scale (estimated from other parameters)\cr
#'  \code{xi}:        \tab MLE of GPD shape\cr
#'  \code{phiu}:      \tab MLE of tail fraction (bulk model or parameterised approach)\cr
#'  \code{se.phiu}:   \tab standard error of MLE of tail fraction\cr
#' }
#' 
#' @note When \code{pvector=NULL} then the initial values are:
#' \itemize{
#'  \item MLE of normal parameters assuming entire population is normal; and
#'  \item threshold 90\% quantile (not relevant for profile likelihood for threshold or fixed threshold approaches);
#'  \item MLE of GPD shape parameter above threshold. 
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
#' @seealso \code{\link[stats:Normal]{dnorm}},
#'  \code{\link[evmix:fgpd]{fgpd}} and \code{\link[evmix:gpd]{gpd}}
#'  
#' @aliases fnormgpdcon lnormgpdcon nlnormgpdcon proflunormgpdcon nlunormgpdcon
#' @family  normgpd normgpdcon gng gngcon fnormgpd fnormgpdcon fgng fgngcon
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
#' # Continuity constraint
#' fit = fnormgpdcon(x)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-4, 4))
#' lines(xx, y)
#' with(fit, lines(xx, dnormgpdcon(xx, nmean, nsd, u, xi), col="red"))
#' abline(v = fit$u, col = "red")
#'   
#' # No continuity constraint
#' fit2 = fnormgpd(x)
#' with(fit2, lines(xx, dnormgpd(xx, nmean, nsd, u, sigmau, xi), col="blue"))
#' abline(v = fit2$u, col = "blue")
#' legend("topleft", c("True Density","No continuity constraint","With continuty constraint"),
#'   col=c("black", "blue", "red"), lty = 1)
#'   
#' # Profile likelihood for initial value of threshold and fixed threshold approach
#' fitu = fnormgpdcon(x, useq = seq(0, 3, length = 20))
#' fitfix = fnormgpdcon(x, useq = seq(0, 3, length = 20), fixedu = TRUE)
#' 
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-4, 4))
#' lines(xx, y)
#' with(fit, lines(xx, dnormgpdcon(xx, nmean, nsd, u, xi), col="red"))
#' abline(v = fit$u, col = "red")
#' with(fitu, lines(xx, dnormgpdcon(xx, nmean, nsd, u, xi), col="purple"))
#' abline(v = fitu$u, col = "purple")
#' with(fitfix, lines(xx, dnormgpdcon(xx, nmean, nsd, u, xi), col="darkgreen"))
#' abline(v = fitfix$u, col = "darkgreen")
#' legend("topleft", c("True Density","Default initial value (90% quantile)",
#'  "Prof. lik. for initial value", "Prof. lik. for fixed threshold"),
#'  col=c("black", "red", "purple", "darkgreen"), lty = 1)
#' }
#'   

# maximum likelihood fitting for normal bulk with GPD for upper tail with single continuity constraint
fnormgpdcon <- function(x, phiu = TRUE, useq = NULL, fixedu = FALSE, pvector = NULL,
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

  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    x = x[is.finite(x)] # ignore missing and infinite cases
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
      pvector[1] = mean(x, trim = 0.2)
      pvector[2] = sd(x)
      pvector[3] = as.vector(quantile(x, 0.9))
      initfgpd = fgpd(x, pvector[3], std.err = FALSE)
      pvector[4] = initfgpd$xi
    }
    
  } else { # profile or fixed
    
    check.nparam(pvector, nparam = np - 1, allownull = TRUE)

    # profile likelihood for threshold or scalar given
    if (length(useq) != 1) {
      
      # remove thresholds with less than 5 excesses
      useq = useq[sapply(useq, FUN = function(u, x) sum(x > u) > 5, x = x)]
      check.param(useq, allowvec = TRUE)
      
      nllhu = sapply(useq, proflunormgpdcon, pvector = pvector, x = x, phiu = phiu,
        method = method, control = control, finitelik = finitelik, ...)
      
      if (all(!is.finite(nllhu))) stop("thresholds are all invalid")
      u = useq[which.min(nllhu)]

    } else {
      u = useq
    }

    if (fixedu) { # threshold fixed
      if (is.null(pvector)) {
        pvector[1] = mean(x, trim = 0.2)
        pvector[2] = sd(x)
        initfgpd = fgpd(x, u, std.err = FALSE)
        pvector[3] = initfgpd$xi
      }
    } else { # threshold as initial value in usual MLE
      if (is.null(pvector)) {
        pvector[1] = mean(x, trim = 0.2)
        pvector[2] = sd(x)
        pvector[3] = u
        initfgpd = fgpd(x, pvector[3], std.err = FALSE)
        pvector[4] = initfgpd$xi
      } else {
        pvector[4] = pvector[3] # shift GPD shape to add in u
        pvector[3] = u
      }
    }
  }

  if (fixedu) { # fixed threshold (separable) likelihood
    nllh = nlunormgpdcon(pvector, u, x, phiu)
    if (is.infinite(nllh)) {
      pvector[3] = 0.1
      nllh = nlunormgpdcon(pvector, u, x, phiu)    
    }
    if (is.infinite(nllh)) stop("initial parameter values are invalid")
  
    fit = optim(par = as.vector(pvector), fn = nlunormgpdcon, u = u, x = x, phiu = phiu,
      finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)    
    
    nmean = fit$par[1]
    nsd = fit$par[2]
    xi = fit$par[3]
    
  } else { # complete (non-separable) likelihood

    nllh = nlnormgpdcon(pvector, x, phiu)
    if (is.infinite(nllh)) {
      pvector[4] = 0.1
      nllh = nlnormgpdcon(pvector, x, phiu)    
    }
    if (is.infinite(nllh)) stop("initial parameter values are invalid")

    fit = optim(par = as.vector(pvector), fn = nlnormgpdcon, x = x, phiu = phiu,
      finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)    
    
    nmean = fit$par[1]
    nsd = fit$par[2]
    u = fit$par[3]
    xi = fit$par[4]
  }
  
  conv = TRUE
  if ((fit$convergence != 0) | any(fit$par == pvector) | (abs(fit$value) >= 1e6)) {
    conv = FALSE
    warning("check convergence")
  }

  pu = pnorm(u, nmean, nsd)
  if (phiu) {
    phiu = 1 - pu
    se.phiu = NA
  } else {
    phiu = mean(x > u, na.rm = TRUE)
    se.phiu = sqrt(phiu * (1 - phiu) / n)
  }
  phib = (1 - phiu) / pu

  du = dnorm(u, nmean, nsd)
  sigmau = phiu / (phib * du)

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
  
  list(call = call, x = as.vector(x), init = as.vector(pvector),
    fixedu = fixedu, useq = useq, nllhuseq = nllhu, optim = fit, conv = conv, cov = invhess,
    mle = fit$par, se = se, rate = phiu, nllh = fit$value, n = n,
    nmean = nmean, nsd = nsd, u = u, sigmau = sigmau, xi = xi, phiu = phiu, se.phiu = se.phiu)
}

#' @export
#' @aliases fnormgpdcon lnormgpdcon nlnormgpdcon proflunormgpdcon nlunormgpdcon
#' @rdname  fnormgpdcon

# log-likelihood function for normal bulk with GPD for upper tail with single continuity constraint
lnormgpdcon <- function(x, nmean = 0, nsd = 1, u = qnorm(0.9, nmean, nsd), 
  xi = 0, phiu = TRUE, log = TRUE) {

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(nmean)
  check.param(nsd)
  check.param(u)
  check.param(xi)
  check.phiu(phiu, allowfalse = TRUE)
  check.logic(log)

  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    x = x[is.finite(x)] # ignore missing and infinite cases
  }

  check.quant(x)
  n = length(x)

  check.inputn(c(length(nmean), length(nsd), length(u), length(xi), length(phiu)), allowscalar = TRUE)

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

  if ((nsd <= 0) | (u <= min(x)) | (u >= max(x))) {
    l = -Inf
  } else {
    if (is.logical(phiu)) {
      pu = pnorm(u, nmean, nsd)
      if (phiu) {
        phiu = 1 - pu
      } else {
        phiu = nu / n
      }
    }
    phib = (1 - phiu) / pu
  
    du = dnorm(u, nmean, nsd)
    sigmau = phiu / (phib * du)

    syu = 1 + xi * (xu - u) / sigmau  
    yb = (xb - nmean) / nsd # used for normal

    
    if ((min(syu) <= 0) | (sigmau <= 0) | (du < .Machine$double.eps) | (phiu <= 0) | (phiu >= 1) | (pu <= 0) | (pu >= 1)) {
      l = -Inf
    } else { 
      l = lgpd(xu, u, sigmau, xi, phiu)
      l = l - nb * log(2 * pi * nsd ^ 2) / 2 - sum(yb ^ 2) / 2 + nb * log(phib)
    }
  }
  
  if (!log) l = exp(l)
  
  l
}

#' @export
#' @aliases fnormgpdcon lnormgpdcon nlnormgpdcon proflunormgpdcon nlunormgpdcon
#' @rdname  fnormgpdcon

# negative log-likelihood function for normal bulk with GPD for upper tail with single continuity constraint
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlnormgpdcon <- function(pvector, x, phiu = TRUE, finitelik = FALSE) {

  np = 4 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.phiu(phiu, allowfalse = TRUE)
  check.logic(finitelik)

  nmean = pvector[1]
  nsd = pvector[2]
  u = pvector[3]
  xi = pvector[4]

  nllh = -lnormgpdcon(x, nmean, nsd, u, xi, phiu) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}

#' @export
#' @aliases fnormgpdcon lnormgpdcon nlnormgpdcon proflunormgpdcon nlunormgpdcon
#' @rdname  fnormgpdcon

# profile negative log-likelihood function for given threshold for
# normal bulk with GPD for upper tail with single continuity constraint
# designed for sapply to loop over vector of thresholds (hence u is first input)
proflunormgpdcon <- function(u, pvector, x, phiu = TRUE, method = "BFGS",
  control = list(maxit = 10000), finitelik = TRUE, ...) {

  np = 4 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np - 1, allownull = TRUE)
  check.param(u)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.phiu(phiu, allowfalse = TRUE)
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
    nllh = nlunormgpdcon(pvector, u, x, phiu)
    
    if (is.infinite(nllh)) pvector = NULL
  }

  if (is.null(pvector)) {
    pvector[1] = mean(x, trim = 0.2)
    pvector[2] = sd(x)
    initfgpd = fgpd(x, u, std.err = FALSE)
    pvector[3] = initfgpd$xi
    nllh = nlunormgpdcon(pvector, u, x, phiu)
  }  

  if (is.infinite(nllh)) {
    pvector[3] = 0.1
    nllh = nlunormgpdcon(pvector, u, x, phiu)    
  }

  # if still invalid then output cleanly
  if (is.infinite(nllh)) {
    warning(paste("initial parameter values for threshold u =", u, "are invalid"))
    fit = list(par = rep(NA, np), value = Inf, counts = 0, convergence = NA, 
      message = "initial values invalid", hessian = rep(NA, np))
  } else {

    fit = optim(par = as.vector(pvector), fn = nlunormgpdcon, u = u, x = x, phiu = phiu,
    finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)
  }
    
  if (finitelik & is.infinite(fit$value)) {
    fit$value = sign(fit$value) * 1e6
  }

  fit$value
}

#' @export
#' @aliases fnormgpdcon lnormgpdcon nlnormgpdcon proflunormgpdcon nlunormgpdcon
#' @rdname  fnormgpdcon

# negative log-likelihood function for normal bulk with GPD for upper tail with single continuity constraint
# (wrapper for likelihood, designed for threshold to be fixed and other parameters optimised)
nlunormgpdcon <- function(pvector, u, x, phiu = TRUE, finitelik = FALSE) {

  np = 4 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np - 1)
  check.param(u)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.phiu(phiu, allowfalse = TRUE)
  check.logic(finitelik)
    
  nmean = pvector[1]
  nsd = pvector[2]
  xi = pvector[3]

  nllh = -lnormgpdcon(x, nmean, nsd, u, xi, phiu) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}
