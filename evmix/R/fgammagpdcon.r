#' @export
#' 
#' @title MLE Fitting of Gamma Bulk and GPD Tail Extreme Value Mixture Model with Single Continuity Constraint
#'
#' @description Maximum likelihood estimation for fitting the extreme value 
#' mixture model with gamma for bulk distribution upto the threshold and conditional
#' GPD above threshold with continuity at threshold. With options for profile likelihood estimation for threshold and
#' fixed threshold approach.
#'
#' @inheritParams fgammagpd
#' 
#' @details The extreme value mixture model with gamma bulk and GPD tail with continuity at threshold is 
#' fitted to the entire dataset using maximum likelihood estimation. The estimated
#' parameters, variance-covariance matrix and their standard errors are automatically
#' output.
#' 
#' See help for \code{\link[evmix:fnormgpd]{fnormgpd}} for details, type \code{help fnormgpd}. 
#' Only the different features are outlined below for brevity.
#' 
#' The GPD \code{sigmau} parameter is now specified as function of other parameters, see 
#' help for \code{\link[evmix:gammagpdcon]{dgammagpdcon}} for details, type \code{help gammagpdcon}.
#' Therefore, \code{sigmau} should not be included in the parameter vector if initial values
#' are provided, making the full parameter vector 
#' (\code{gshape}, \code{gscale}, \code{u}, \code{xi}) if threshold is also estimated and
#' (\code{gshape}, \code{gscale}, \code{xi}) for profile likelihood or fixed threshold approach.
#' 
#' Non-positive data are ignored as likelihood is infinite, except for \code{gshape=1}.
#' 
#' @return Log-likelihood is given by \code{\link[evmix:fgammagpdcon]{lgammagpdcon}} and it's
#'   wrappers for negative log-likelihood from \code{\link[evmix:fgammagpdcon]{nlgammagpdcon}}
#'   and \code{\link[evmix:fgammagpdcon]{nlugammagpdcon}}. Profile likelihood for single
#'   threshold given by \code{\link[evmix:fgammagpdcon]{proflugammagpdcon}}. Fitting function
#'   \code{\link[evmix:fgammagpdcon]{fgammagpdcon}} returns a simple list with the
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
#'  \code{gshape}:    \tab MLE of gamma shape\cr
#'  \code{gscale}:    \tab MLE of gamma scale\cr
#'  \code{u}:         \tab threshold (fixed or MLE)\cr
#'  \code{sigmau}:    \tab MLE of GPD scale (estimated from other parameters)\cr
#'  \code{xi}:        \tab MLE of GPD shape\cr
#'  \code{phiu}:      \tab MLE of tail fraction (bulk model or parameterised approach)\cr
#'  \code{se.phiu}:   \tab standard error of MLE of tail fraction\cr
#' }
#' 
#' @note When \code{pvector=NULL} then the initial values are:
#' \itemize{
#'  \item approximation of MLE of gamma parameters assuming entire population is gamma; and
#'  \item threshold 90\% quantile (not relevant for profile likelihood for threshold or fixed threshold approaches);
#'  \item MLE of GPD shape parameter above threshold. 
#' }
#' 
#' @references
#' \url{http://www.math.canterbury.ac.nz/~c.scarrott/evmix}
#' 
#' \url{http://en.wikipedia.org/wiki/Gamma_distribution}
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
#' @seealso \code{\link[stats:GammaDist]{dgamma}},
#'  \code{\link[evmix:fgpd]{fgpd}} and \code{\link[evmix:gpd]{gpd}}
#'  
#' @aliases fgammagpdcon lgammagpdcon nlgammagpdcon proflugammagpdcon nlugammagpdcon
#' @family  mgamma fmgamma
#'          gammagpd gammagpdcon fgammagpd fgammagpdcon normgpd fnormgpd
#'          mgammagpd mgammagpdcon fmgammagpd fmgammagpdcon 
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' par(mfrow = c(2, 1))
#' 
#' x = rgamma(1000, shape = 2)
#' xx = seq(-0.1, 8, 0.01)
#' y = dgamma(xx, shape = 2)
#' 
#' # Continuity constraint
#' fit = fgammagpdcon(x)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-0.1, 8))
#' lines(xx, y)
#' with(fit, lines(xx, dgammagpdcon(xx, gshape, gscale, u, xi), col="red"))
#' abline(v = fit$u, col = "red")
#'   
#' # No continuity constraint
#' fit2 = fgammagpd(x, phiu = FALSE)
#' with(fit2, lines(xx, dgammagpd(xx, gshape, gscale, u, sigmau, xi, phiu), col="blue"))
#' abline(v = fit2$u, col = "blue")
#' legend("topright", c("True Density","No continuity constraint","With continuty constraint"),
#'   col=c("black", "blue", "red"), lty = 1)
#'   
#' # Profile likelihood for initial value of threshold and fixed threshold approach
#' fitu = fgammagpdcon(x, useq = seq(1, 5, length = 20))
#' fitfix = fgammagpdcon(x, useq = seq(1, 5, length = 20), fixedu = TRUE)
#' 
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-0.1, 8))
#' lines(xx, y)
#' with(fit, lines(xx, dgammagpdcon(xx, gshape, gscale, u, xi), col="red"))
#' abline(v = fit$u, col = "red")
#' with(fitu, lines(xx, dgammagpdcon(xx, gshape, gscale, u, xi), col="purple"))
#' abline(v = fitu$u, col = "purple")
#' with(fitfix, lines(xx, dgammagpdcon(xx, gshape, gscale, u, xi), col="darkgreen"))
#' abline(v = fitfix$u, col = "darkgreen")
#' legend("topright", c("True Density","Default initial value (90% quantile)",
#'  "Prof. lik. for initial value", "Prof. lik. for fixed threshold"),
#'  col=c("black", "red", "purple", "darkgreen"), lty = 1)
#' }
#'   

# maximum likelihood fitting for gamma bulk with GPD for upper tail with continuity at threshold
fgammagpdcon <- function(x, phiu = TRUE, useq = NULL, fixedu = FALSE, pvector = NULL,
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
      s = log(mean(x)) - mean(log(x))
      k = (3 - s + sqrt( (s-3) ^ 2 + 24 * s)) / 12 / s
      pvector[1] = k
      pvector[2] = mean(x) / k
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
      check.posparam(useq, allowvec = TRUE)
      
      nllhu = sapply(useq, proflugammagpdcon, pvector = pvector, x = x, phiu = phiu,
        method = method, control = control, finitelik = finitelik, ...)
      
      if (all(!is.finite(nllhu))) stop("thresholds are all invalid")
      u = useq[which.min(nllhu)]

    } else {
      u = useq
    }

    if (fixedu) { # threshold fixed
      if (is.null(pvector)) {
        s = log(mean(x)) - mean(log(x))
        k = (3 - s + sqrt( (s-3) ^ 2 + 24 * s)) / 12 / s
        pvector[1] = k
        pvector[2] = mean(x) / k
        initfgpd = fgpd(x, u, std.err = FALSE)
        pvector[3] = initfgpd$xi
      }
    } else { # threshold as initial value in usual MLE
      if (is.null(pvector)) {
        s = log(mean(x)) - mean(log(x))
        k = (3 - s + sqrt( (s-3) ^ 2 + 24 * s)) / 12 / s
        pvector[1] = k
        pvector[2] = mean(x) / k
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
    nllh = nlugammagpdcon(pvector, u, x, phiu)
    if (is.infinite(nllh)) {
      pvector[3] = 0.1
      nllh = nlugammagpdcon(pvector, u, x, phiu)    
    }
    if (is.infinite(nllh)) stop("initial parameter values are invalid")
  
    fit = optim(par = as.vector(pvector), fn = nlugammagpdcon, u = u, x = x, phiu = phiu,
      finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)    
    
    gshape = fit$par[1]
    gscale = fit$par[2]
    xi = fit$par[3]
    
  } else { # complete (non-separable) likelihood
    
    nllh = nlgammagpdcon(pvector, x, phiu)
    if (is.infinite(nllh)) {
      pvector[4] = 0.1
      nllh = nlgammagpdcon(pvector, x, phiu)    
    }
    if (is.infinite(nllh)) stop("initial parameter values are invalid")
  
    fit = optim(par = as.vector(pvector), fn = nlgammagpdcon, x = x, phiu = phiu,
      finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)    
    
    gshape = fit$par[1]
    gscale = fit$par[2]
    u = fit$par[3]
    xi = fit$par[4]
  }
  
  conv = TRUE
  if ((fit$convergence != 0) | any(fit$par == pvector) | (abs(fit$value) >= 1e6)) {
    conv = FALSE
    warning("check convergence")
  }

  pu = pgamma(u, gshape, scale = gscale)
  if (phiu) {
    phiu = 1 - pu
    se.phiu = NA
  } else {
    phiu = mean(x > u, na.rm = TRUE)
    se.phiu = sqrt(phiu * (1 - phiu) / n)
  }
  phib = (1 - phiu) / pu

  du = dgamma(u, gshape, scale = gscale)
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

  list(call = call, x = as.vector(x), 
    init = as.vector(pvector), fixedu = fixedu, useq = useq, nllhuseq = nllhu,
    optim = fit, conv = conv, cov = invhess, mle = fit$par, se = se, rate = phiu,
    nllh = fit$value, n = n,
    gshape = gshape, gscale = gscale, u = u, sigmau = sigmau, xi = xi, phiu = phiu, se.phiu = se.phiu)
}

#' @export
#' @aliases fgammagpdcon lgammagpdcon nlgammagpdcon proflugammagpdcon nlugammagpdcon
#' @rdname  fgammagpdcon

# log-likelihood function for gamma bulk with GPD for upper tail with continuity at threshold
lgammagpdcon <- function(x, gshape = 1, gscale = 1, u = qgamma(0.9, gshape, 1/gscale),
  xi = 0, phiu = TRUE, log = TRUE) {

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(gshape)
  check.param(gscale)
  check.param(u)
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

  check.inputn(c(length(gshape), length(gscale), length(u), length(xi), length(phiu)), allowscalar = TRUE)

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

  if ((gscale <= 0) | (gshape <= 0) | (u <= 0) | (u <= min(x)) | (u >= max(x))) {
    l = -Inf
  } else {
    if (is.logical(phiu)) {
      pu = pgamma(u, gshape, scale = gscale)
      if (phiu) {
        phiu = 1 - pu
      } else {
        phiu = nu / n
      }
    }
    phib = (1 - phiu) / pu
  
    du = dgamma(u, gshape, scale = gscale)
    sigmau = phiu / (phib * du)
    
    syu = 1 + xi * (xu - u) / sigmau  
  
    if ((min(syu) <= 0) | (sigmau <= 0) | (du < .Machine$double.eps) | (phiu <= 0) | (phiu >= 1) | (pu <= 0) | (pu >= 1)) {
      l = -Inf
    } else { 
      l = lgpd(xu, u, sigmau, xi, phiu)
      l = l + (gshape - 1) * sum(log(xb)) - sum(xb) / gscale - nb * gshape * log(gscale) - nb * lgamma(gshape)  + nb * log(phib)
    }
  }
  
  if (!log) l = exp(l)
  
  l
}

#' @export
#' @aliases fgammagpdcon lgammagpdcon nlgammagpdcon proflugammagpdcon nlugammagpdcon
#' @rdname  fgammagpdcon

# negative log-likelihood function for gamma bulk with GPD for upper tail with continuity at threshold
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlgammagpdcon <- function(pvector, x, phiu = TRUE, finitelik = FALSE) {

  np = 4 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.phiu(phiu, allowfalse = TRUE)
  check.logic(finitelik)

  gshape = pvector[1]
  gscale = pvector[2]
  u = pvector[3]
  xi = pvector[4]

  nllh = -lgammagpdcon(x, gshape, gscale, u, xi, phiu) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}

#' @export
#' @aliases fgammagpdcon lgammagpdcon nlgammagpdcon proflugammagpdcon nlugammagpdcon
#' @rdname  fgammagpdcon

# profile negative log-likelihood function for given threshold for
# gamma bulk with GPD for upper tail with continuity at threshold
# designed for sapply to loop over vector of thresholds (hence u is first input)
proflugammagpdcon <- function(u, pvector, x, phiu = TRUE, method = "BFGS",
  control = list(maxit = 10000), finitelik = TRUE, ...) {

  np = 4 # maximum number of parameters

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
    nllh = nlugammagpdcon(pvector, u, x, phiu)
    
    if (is.infinite(nllh)) pvector = NULL
  }

  if (is.null(pvector)) {
    s = log(mean(x)) - mean(log(x))
    k = (3 - s + sqrt( (s-3) ^ 2 + 24 * s)) / 12 / s
    pvector[1] = k
    pvector[2] = mean(x) / k
    initfgpd = fgpd(x, u, std.err = FALSE)
    pvector[3] = initfgpd$xi
    nllh = nlugammagpdcon(pvector, u, x, phiu)
  }  

  if (is.infinite(nllh)) {
    pvector[3] = 0.1
    nllh = nlugammagpdcon(pvector, u, x, phiu)    
  }

  # if still invalid then output cleanly
  if (is.infinite(nllh)) {
    warning(paste("initial parameter values for threshold u =", u, "are invalid"))
    fit = list(par = rep(NA, np), value = Inf, counts = 0, convergence = NA, 
      message = "initial values invalid", hessian = rep(NA, np))
  } else {

    fit = optim(par = as.vector(pvector), fn = nlugammagpdcon, u = u, x = x, phiu = phiu,
    finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)
  }
    
  if (finitelik & is.infinite(fit$value)) {
    fit$value = sign(fit$value) * 1e6
  }

  fit$value
}

#' @export
#' @aliases fgammagpdcon lgammagpdcon nlgammagpdcon proflugammagpdcon nlugammagpdcon
#' @rdname  fgammagpdcon

# negative log-likelihood function for gamma bulk with GPD for upper tail with continuity at threshold
# (wrapper for likelihood, designed for threshold to be fixed and other parameters optimised)
nlugammagpdcon <- function(pvector, u, x, phiu = TRUE, finitelik = FALSE) {

  np = 4 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np - 1)
  check.posparam(u)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.phiu(phiu, allowfalse = TRUE)
  check.logic(finitelik)
    
  gshape = pvector[1]
  gscale = pvector[2]
  xi = pvector[3]

  nllh = -lgammagpdcon(x, gshape, gscale, u, xi, phiu) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}
