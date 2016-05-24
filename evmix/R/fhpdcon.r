#' @export
#' 
#' @title MLE Fitting of Hybrid Pareto Extreme Value Mixture Model with Single Continuity Constraint
#'
#' @description Maximum likelihood estimation for fitting the Hybrid Pareto extreme
#' value mixture model, with only continuity at threshold and not necessarily
#' continuous in first derivative. With options for profile likelihood estimation for
#' threshold and fixed threshold approach.
#'
#' @inheritParams fnormgpd
#' @inheritParams dnormgpd
#' @inheritParams fgpd
#' 
#' @details The hybrid Pareto model is fitted to the entire dataset using maximum
#' likelihood estimation, with only continuity at threshold and not necessarily
#' continuous in first derivative. The estimated parameters, variance-covariance matrix
#' and their standard errors are automatically output.
#' 
#' Note that the key difference between this model (\code{hpdcon}) and the 
#' normal with GPD tail and continuity at threshold (\code{normgpdcon}) is that the
#' latter includes the rescaling of the conditional GPD component
#' by the tail fraction to make it an unconditional tail model. However, for the hybrid
#' Pareto with single continuity constraint use the GPD in it's conditional form with no
#' differential scaling compared to the bulk model.
#' 
#' See help for \code{\link[evmix:fnormgpd]{fnormgpd}} for details, type \code{help fnormgpd}. Only
#' the different features are outlined below for brevity.
#' 
#' The profile likelihood and fixed threshold approach functionality are implemented for this
#' version of the hybrid Pareto as it includes the threshold as a parameter. Whereas the usual
#' hybrid Pareto does not naturally have a threshold parameter.
#' 
#' The GPD \code{sigmau} parameter is now specified as function of other parameters, see 
#' help for \code{\link[evmix:hpdcon]{dhpdcon}} for details, type \code{help hpdcon}.
#' Therefore, \code{sigmau} should not be included in the parameter vector if initial values
#' are provided, making the full parameter vector 
#' (\code{nmean}, \code{nsd}, \code{u}, \code{xi}) if threshold is also estimated and
#' (\code{nmean}, \code{nsd}, \code{xi}) for profile likelihood or fixed threshold approach.
#' 
#' @return \code{\link[evmix:fhpdcon]{lhpdcon}}, \code{\link[evmix:fhpdcon]{nlhpdcon}},
#' and \code{\link[evmix:fhpdcon]{nluhpdcon}} give the log-likelihood,
#' negative log-likelihood and profile likelihood for threshold. Profile likelihood
#' for single threshold is given by \code{\link[evmix:fhpdcon]{profluhpdcon}}.
#' \code{\link[evmix:fhpdcon]{fhpdcon}} returns a simple list with the following elements
#'
#' \tabular{ll}{
#'  \code{call}:    \tab \code{optim} call\cr
#'  \code{x}:       \tab data vector \code{x}\cr
#'  \code{init}:    \tab \code{pvector}\cr
#'  \code{fixedu}:  \tab fixed threshold, logical\cr
#'  \code{useq}:    \tab threshold vector for profile likelihood or scalar for fixed threshold\cr
#'  \code{nllhuseq}:  \tab profile negative log-likelihood at each threshold in useq\cr
#'  \code{optim}:   \tab complete \code{optim} output\cr
#'  \code{mle}:     \tab vector of MLE of parameters\cr
#'  \code{cov}:     \tab variance-covariance matrix of MLE of parameters\cr
#'  \code{se}:      \tab vector of standard errors of MLE of parameters\cr
#'  \code{rate}:    \tab \code{phiu} to be consistent with \code{\link[evd:fpot]{evd}}\cr
#'  \code{nllh}:    \tab minimum negative log-likelihood\cr
#'  \code{n}:       \tab total sample size\cr
#'  \code{nmean}:   \tab MLE of normal mean\cr
#'  \code{nsd}:     \tab MLE of normal standard deviation\cr
#'  \code{u}:       \tab threshold (fixed or MLE)\cr
#'  \code{sigmau}:  \tab MLE of GPD scale (estimated from other parameters)\cr
#'  \code{xi}:      \tab MLE of GPD shape\cr
#'  \code{phiu}:    \tab MLE of tail fraction \code{1/(1+pnorm(u,nmean,nsd))}\cr
#' }
#' 
#' @note When \code{pvector=NULL} then the initial values are:
#' \itemize{
#'  \item threshold 90\% quantile (not relevant for profile likelihood for threshold or fixed threshold approaches);
#'  \item MLE of normal parameters assuming entire population is normal; and
#'  \item MLE of GPD parameters above threshold. 
#' }
#' Avoid setting the starting value for the shape parameter to
#' \code{xi=0} as depending on the optimisation method it may be get stuck.
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
#' Carreau, J. and Y. Bengio (2008). A hybrid Pareto model for asymmetric fat-tailed data:
#' the univariate case. Extremes 12 (1), 53-76.
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @section Acknowledgments: See Acknowledgments in
#'   \code{\link[evmix:fnormgpd]{fnormgpd}}, type \code{help fnormgpd}.
#' 
#' @seealso \code{\link[stats:Normal]{dnorm}},
#'  \code{\link[evmix:fgpd]{fgpd}} and \code{\link[evmix:gpd]{gpd}}
#'  
#' The \code{\link[condmixt:condmixt-package]{condmixt}} package written by one of the
#' original authors of the hybrid Pareto model (Carreau and Bengio, 2008) also has 
#' similar functions for the likelihood of the hybrid Pareto 
#' \code{\link[condmixt:hpareto.negloglike]{hpareto.negloglike}} and
#' fitting \code{\link[condmixt:hpareto.negloglike]{hpareto.fit}}.
#' 
#' @aliases fhpdcon lhpdcon nlhpdcon profluhpdcon nluhpdcon
#' @family  hpd hpdcon normgpd normgpdcon gng gngcon
#'          fhpd fhpdcon fnormgpd fnormgpdcon fgng fgngcon
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
#' # Hybrid Pareto provides reasonable fit for some asymmetric heavy upper tailed distributions
#' # but not for cases such as the normal distribution
#' 
#' # Continuity constraint
#' fit = fhpdcon(x)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-4, 4))
#' lines(xx, y)
#' with(fit, lines(xx, dhpdcon(xx, nmean, nsd, u, xi), col="red"))
#' abline(v = fit$u, col = "red")
#'   
#' # No continuity constraint
#' fit2 = fhpd(x)
#' with(fit2, lines(xx, dhpd(xx, nmean, nsd, xi), col="blue"))
#' abline(v = fit2$u, col = "blue")
#' legend("topleft", c("True Density","No continuity constraint","With continuty constraint"),
#'   col=c("black", "blue", "red"), lty = 1)
#'   
#' # Profile likelihood for initial value of threshold and fixed threshold approach
#' fitu = fhpdcon(x, useq = seq(-2, 2, length = 20))
#' fitfix = fhpdcon(x, useq = seq(-2, 2, length = 20), fixedu = TRUE)
#' 
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-4, 4))
#' lines(xx, y)
#' with(fit, lines(xx, dhpdcon(xx, nmean, nsd, u, xi), col="red"))
#' abline(v = fit$u, col = "red")
#' with(fitu, lines(xx, dhpdcon(xx, nmean, nsd, u, xi), col="purple"))
#' abline(v = fitu$u, col = "purple")
#' with(fitfix, lines(xx, dhpdcon(xx, nmean, nsd, u, xi), col="darkgreen"))
#' abline(v = fitfix$u, col = "darkgreen")
#' legend("topleft", c("True Density","Default initial value (90% quantile)",
#'  "Prof. lik. for initial value", "Prof. lik. for fixed threshold"),
#'  col=c("black", "red", "purple", "darkgreen"), lty = 1)
#'   
#' # Notice that if tail fraction is included a better fit is obtained
#' fittailfrac = fnormgpdcon(x)
#' 
#' par(mfrow = c(1, 1))
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-4, 4))
#' lines(xx, y)
#' with(fit, lines(xx, dhpdcon(xx, nmean, nsd, u, xi), col="red"))
#' abline(v = fit$u, col = "red")
#' with(fittailfrac, lines(xx, dnormgpdcon(xx, nmean, nsd, u, xi), col="blue"))
#' abline(v = fittailfrac$u)
#' legend("topright", c("Standard Normal", "Hybrid Pareto Continuous", "Normal+GPD Continuous"),
#'   col=c("black", "red", "blue"), lty = 1)
#' }
#' 

# maximum likelihood fitting for hybrid Pareto with single continuity constraint
fhpdcon <- function(x, useq = NULL, fixedu = FALSE, pvector = NULL,
  std.err = TRUE, method = "BFGS", control = list(maxit = 10000), finitelik = TRUE, ...) {

  call <- match.call()
    
  np = 4 # maximum number of parameters

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
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
      
      nllhu = sapply(useq, profluhpdcon, pvector = pvector, x = x,
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
    nllh = nluhpdcon(pvector, u, x)
    if (is.infinite(nllh)) {
      pvector[3] = 0.1
      nllh = nluhpdcon(pvector, u, x)    
    }
    if (is.infinite(nllh)) stop("initial parameter values are invalid")
  
    fit = optim(par = as.vector(pvector), fn = nluhpdcon, u = u, x = x,
      finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)    
    
    nmean = fit$par[1]
    nsd = fit$par[2]
    xi = fit$par[3]
    
  } else { # complete (non-separable) likelihood

    nllh = nlhpdcon(pvector, x)
    if (is.infinite(nllh)) {
      pvector[4] = 0.1
      nllh = nlhpdcon(pvector, x)    
    }
    if (is.infinite(nllh)) stop("initial parameter values are invalid")

    fit = optim(par = as.vector(pvector), fn = nlhpdcon, x = x,
      finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)    
    
    nmean = fit$par[1]
    nsd = fit$par[2]
    u = fit$par[3]
    xi = fit$par[4]
  }
  
  du = dnorm(u, nmean, nsd)
  sigmau = 1/du

  r = 1 + pnorm(u, nmean, nsd)
  phiu = 1/r # same normalisation for GPD and normal
  
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

  list(call = call, x = as.vector(x), 
    init = as.vector(pvector), fixedu = fixedu, useq = useq, nllhuseq = nllhu,
    optim = fit, conv = conv, cov = invhess, mle = fit$par, se = se, rate = phiu,
    nllh = fit$value, n = n,
    nmean = nmean, nsd = nsd, u = u, sigmau = sigmau, xi = xi, phiu = phiu)
}

#' @export
#' @aliases fhpdcon lhpdcon nlhpdcon profluhpdcon nluhpdcon
#' @rdname  fhpdcon

# log-likelihood function for hybrid Pareto with single continuity constraint
lhpdcon <- function(x, nmean = 0, nsd = 1, u = qnorm(0.9, nmean, nsd), 
  xi = 0, log = TRUE) {

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(nmean)
  check.param(nsd)
  check.param(u)
  check.param(xi)
  check.logic(log)

  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    x = x[is.finite(x)] # ignore missing and infinite cases
  }

  check.quant(x)
  n = length(x)

  check.inputn(c(length(nmean), length(nsd), length(u), length(xi)), allowscalar = TRUE)

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
    du = dnorm(u, nmean, nsd)
    sigmau = 1/du
  
    r = 1 + pnorm(u, nmean, nsd)

    syu = 1 + xi * (xu - u) / sigmau  
    yb = (xb - nmean) / nsd # used for normal
  
    if ((min(syu) <= 0) | (sigmau <= 0) | (du < .Machine$double.eps)) {
      l = -Inf
    } else { 
      l = lgpd(xu, u, sigmau, xi) # phiu disappears
      l = l - nb * log(2 * pi * nsd ^ 2) / 2 - sum(yb ^ 2) / 2 # phib disappears
      l = l - n * log(r) # divide by normalisation constant
    }
  }
  
  if (!log) l = exp(l)
  
  l
}

#' @export
#' @aliases fhpdcon lhpdcon nlhpdcon profluhpdcon nluhpdcon
#' @rdname  fhpdcon

# negative log-likelihood function for hybrid Pareto with single continuity constraint
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlhpdcon <- function(pvector, x, finitelik = FALSE) {

  np = 4 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.logic(finitelik)

  nmean = pvector[1]
  nsd = pvector[2]
  u = pvector[3]
  xi = pvector[4]

  nllh = -lhpdcon(x, nmean, nsd, u, xi) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}

#' @export
#' @aliases fhpdcon lhpdcon nlhpdcon profluhpdcon nluhpdcon
#' @rdname  fhpdcon

# profile negative log-likelihood function for given threshold for
# hybrid Pareto with single continuity constraint
# designed for sapply to loop over vector of thresholds (hence u is first input)
profluhpdcon <- function(u, pvector, x,
  method = "BFGS", control = list(maxit = 10000), finitelik = TRUE, ...) {

  np = 4 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np - 1, allownull = TRUE)
  check.param(u)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.optim(method)
  check.control(control)
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
    nllh = nluhpdcon(pvector, u, x)
    
    if (is.infinite(nllh)) pvector = NULL
  }

  if (is.null(pvector)) {
    pvector[1] = mean(x, trim = 0.2)
    pvector[2] = sd(x)
    initfgpd = fgpd(x, u, std.err = FALSE)
    pvector[3] = initfgpd$xi
    nllh = nluhpdcon(pvector, u, x)
  }  

  if (is.infinite(nllh)) {
    pvector[3] = 0.1
    nllh = nluhpdcon(pvector, u, x)    
  }

  # if still invalid then output cleanly
  if (is.infinite(nllh)) {
    warning(paste("initial parameter values for threshold u =", u, "are invalid"))
    fit = list(par = rep(NA, np), value = Inf, counts = 0, convergence = NA, 
      message = "initial values invalid", hessian = rep(NA, np))
  } else {

    fit = optim(par = as.vector(pvector), fn = nluhpdcon, u = u, x = x,
    finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)
  }
    
  if (finitelik & is.infinite(fit$value)) {
    fit$value = sign(fit$value) * 1e6
  }

  fit$value
}

#' @export
#' @aliases fhpdcon lhpdcon nlhpdcon profluhpdcon nluhpdcon
#' @rdname  fhpdcon

# negative log-likelihood function for hybrid Pareto with single continuity constraint
# (wrapper for likelihood, designed for threshold to be fixed and other parameters optimised)
nluhpdcon <- function(pvector, u, x, finitelik = FALSE) {

  np = 4 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np - 1)
  check.param(u)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.logic(finitelik)
    
  nmean = pvector[1]
  nsd = pvector[2]
  xi = pvector[3]

  nllh = -lhpdcon(x, nmean, nsd, u, xi) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}
