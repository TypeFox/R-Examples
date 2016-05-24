#' @export
#' 
#' @title MLE Fitting of Normal Bulk and GPD Tail Extreme Value Mixture Model
#'
#' @description Maximum likelihood estimation for fitting the extreme value 
#' mixture model with normal for bulk distribution upto the threshold and conditional
#' GPD above threshold. With options for profile likelihood estimation for threshold and
#' fixed threshold approach.
#'
#' @param nmean   scalar normal mean
#' @param nsd     scalar normal standard deviation (positive)
#' @param u       scalar threshold value
#' @param pvector vector of initial values of parameters or \code{NULL} for default
#'                values, see below
#' @param phiu    probability of being above threshold \eqn{(0, 1)} or logical, see Details in 
#'                help for \code{\link[evmix:fnormgpd]{fnormgpd}}
#' @param useq    vector of thresholds (or scalar) to be considered in profile likelihood or
#'                \code{NULL} for no profile likelihood
#' @param fixedu  logical, should threshold be fixed (at either scalar value in \code{useq},
#'                or estimated from maximum of profile likelihood evaluated at
#'                sequence of thresholds in \code{useq})
#' @inheritParams fgpd
#' 
#' @details The extreme value mixture model with normal bulk and GPD tail is 
#' fitted to the entire dataset using maximum likelihood estimation. The estimated
#' parameters, variance-covariance matrix and their standard errors are automatically
#' output.
#' 
#' The optimisation of the likelihood for these mixture models can be very sensitive to
#' the initial parameter vector (particularly the threshold), as often there are numerous
#' local modes where multiple thresholds give similar fits. This is an inherent feature
#' of such models. Options are provided by the arguments \code{pvector},
#' \code{useq} and \code{fixedu} to implement various commonly used likelihood inference
#' approaches for such models:
#' \enumerate{
#'  \item (default) \code{pvector=NULL}, \code{useq=NULL} and \code{fixedu=FALSE} 
#'    - to set initial value for threshold at 90\% quantile along with usual defaults for
#'    other parameters as defined in Notes below. Standard likelihood optimisation is used;
#'  \item \code{pvector=c(nmean, nsd, u, sigmau, xi)} - where initial values of all
#'    5 parameters are manually set. Standard likelihood optimisation is used;
#'  \item \code{useq} as vector - to specify a sequence of thresholds at which to evaluate
#'    profile likelihood and extract threshold which gives maximum profile likelihood; or
#'  \item \code{useq} as scalar - to specify a single value for threshold to be considered.
#' }
#' In options (3) and (4) the threshold can be treated as: 
#' \itemize{
#'  \item initial value for maximum likelihood estimation when \code{fixedu=FALSE}, using
#'    either profile likelihood estimate (3) or pre-chosen threshold (4); or
#'  \item a fixed threshold with MLE for other parameters when \code{fixedu=TRUE}, using
#'    either profile likelihood estimate (3) or pre-chosen threshold (4).
#' }
#' The latter approach can be used to implement the traditional fixed threshold modelling
#' approach with threshold pre-chosen using, for example, graphical diagnostics. Further,
#' in either such case (3) or (4) the \code{pvector} could be:
#' \itemize{
#'  \item \code{NULL} for usual defaults for other four parameters, defined in Notes below; or
#'  \item vector of initial values for remaining 4 parameters 
#'    (\code{nmean}, \code{nsd}, \code{sigmau}, \code{xi}).
#' }
#' If the threshold is treated as fixed, then the likelihood is separable between the bulk
#' and tail components. However, in practice we have found black-box optimisation of the
#' combined likelihood works sufficiently well, so is used herein.
#' 
#' The following functions are provided:
#' \itemize{
#'  \item \code{\link[evmix:fnormgpd]{fnormgpd}} - maximum likelihood fitting with all the above options;
#'  \item \code{\link[evmix:fnormgpd]{lnormgpd}} - log-likelihood;
#'  \item \code{\link[evmix:fnormgpd]{nlnormgpd}} - negative log-likelihood;
#'  \item \code{\link[evmix:fnormgpd]{proflunormgpd}} - profile likelihood for given threshold; and
#'  \item \code{\link[evmix:fnormgpd]{nlunormgpd}} - negative log-likelihood (threshold specified separately).
#' }
#' The log-likelihood functions are provided for wider usage, e.g. constructing
#' profile likelihood functions.
#' 
#' Defaults values for the parameter vector \code{pvector} are given in the fitting 
#' \code{\link[evmix:fnormgpd]{fnormgpd}} and profile likelihood functions
#' \code{\link[evmix:fnormgpd]{proflunormgpd}}. The parameter vector \code{pvector}
#' must be specified in the negative log-likelihood functions 
#' \code{\link[evmix:fnormgpd]{nlnormgpd}} and \code{\link[evmix:fnormgpd]{nlunormgpd}}. 
#' The threshold \code{u} must also be specified in the profile likelihood function
#' \code{\link[evmix:fnormgpd]{proflunormgpd}} and \code{\link[evmix:fnormgpd]{nlunormgpd}}.
#' 
#' Log-likelihood calculations are carried out in \code{\link[evmix:fnormgpd]{lnormgpd}},
#' which takes parameters as inputs in the same form as distribution functions. The negative
#' log-likelihood functions \code{\link[evmix:fnormgpd]{nlnormgpd}} and
#' \code{\link[evmix:fnormgpd]{nlunormgpd}} are wrappers for likelihood function
#' \code{\link[evmix:fnormgpd]{lnormgpd}} designed towards optimisation, 
#' i.e. \code{\link[evmix:fnormgpd]{nlnormgpd}} has vector of all 5 parameters as
#' first input and \code{\link[evmix:fnormgpd]{nlunormgpd}} has threshold as second input
#' and vector of remaining 4 parameters as first input. The profile likelihood
#' function \code{\link[evmix:fnormgpd]{proflunormgpd}} has threshold \code{u} as the first
#' input, to permit use of \code{\link[base:lapply]{sapply}} function to evaluate profile
#' likelihood over vector of potential thresholds. 
#' 
#' The tail fraction \code{phiu} is treated separately to the other parameters, 
#' to allow for all it's representations. In the fitting 
#' \code{\link[evmix:fnormgpd]{fnormgpd}} and profile likelihood function
#' \code{\link[evmix:fnormgpd]{proflunormgpd}} it is logical:
#' \itemize{
#'  \item default value \code{phiu=TRUE} - tail fraction specified by 
#'    normal survivor function \code{phiu = 1 - pnorm(u, nmean, nsd)} and standard error is
#'    output as \code{NA}; and
#'  \item \code{phiu=FALSE} - treated as extra parameter estimated using the MLE which is
#'    the sample proportion above the threshold and standard error is output.
#' }
#' In the likelihood functions \code{\link[evmix:fnormgpd]{lnormgpd}},
#' \code{\link[evmix:fnormgpd]{nlnormgpd}} and \code{\link[evmix:fnormgpd]{nlunormgpd}} 
#' it can be logical or numeric:
#' \itemize{
#'  \item logical - same as for fitting functions with default value \code{phiu=TRUE}.
#'  \item numeric - any value over range \eqn{(0, 1)}. Notice that the tail
#'    fraction probability cannot be 0 or 1 otherwise there would be no
#'    contribution from either tail or bulk components respectively.
#' }
#' 
#' Missing values (\code{NA} and \code{NaN}) are assumed to be invalid data so are ignored,
#' which is inconsistent with the \code{\link[evd:fpot]{evd}} library which assumes the 
#' missing values are below the threshold.
#' 
#' The function \code{\link[evmix:fnormgpd]{lnormgpd}} carries out the calculations
#' for the log-likelihood directly, which can be exponentiated to give actual
#' likelihood using (\code{log=FALSE}).
#' 
#' The default optimisation algorithm is "BFGS", which requires a finite negative 
#' log-likelihood function evaluation \code{finitelik=TRUE}. For invalid 
#' parameters, a zero likelihood is replaced with \code{exp(-1e6)}. The "BFGS" 
#' optimisation algorithms require finite values for likelihood, so any user 
#' input for \code{finitelik} will be overridden and set to \code{finitelik=TRUE} 
#' if either of these optimisation methods is chosen.
#' 
#' It will display a warning for non-zero convergence result comes from 
#' \code{\link[stats:optim]{optim}} function call or for common indicators of lack
#' of convergence (e.g. any estimated parameters same as initial values).
#' 
#' If the hessian is of reduced rank then the variance covariance (from inverse hessian)
#' and standard error of parameters cannot be calculated, then by default 
#' \code{std.err=TRUE} and the function will stop. If you want the parameter estimates
#' even if the hessian is of reduced rank (e.g. in a simulation study) then
#' set \code{std.err=FALSE}. 
#' 
#' @return Log-likelihood is given by \code{\link[evmix:fnormgpd]{lnormgpd}} and it's
#'   wrappers for negative log-likelihood from \code{\link[evmix:fnormgpd]{nlnormgpd}}
#'   and \code{\link[evmix:fnormgpd]{nlunormgpd}}. Profile likelihood for single
#'   threshold given by \code{\link[evmix:fnormgpd]{proflunormgpd}}. Fitting function
#'   \code{\link[evmix:fnormgpd]{fnormgpd}} returns a simple list with the
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
#'  \code{sigmau}:    \tab MLE of GPD scale\cr
#'  \code{xi}:        \tab MLE of GPD shape\cr
#'  \code{phiu}:      \tab MLE of tail fraction (bulk model or parameterised approach)\cr
#'  \code{se.phiu}:   \tab standard error of MLE of tail fraction\cr
#' }
#' 
#' The output list has some duplicate entries and repeats some of the inputs to both 
#' provide similar items to those from \code{\link[evd:fpot]{fpot}} and increase usability.
#'  
#' @note Unlike most of the distribution functions for the extreme value mixture models,
#' the MLE fitting only permits single scalar values for each parameter and 
#' \code{phiu}.
#' 
#' When \code{pvector=NULL} then the initial values are:
#' \itemize{
#'  \item MLE of normal parameters assuming entire population is normal; and
#'  \item threshold 90\% quantile (not relevant for profile likelihood or fixed threshold approaches);
#'  \item MLE of GPD parameters above threshold. 
#' }
#' Avoid setting the starting value for the shape parameter to
#' \code{xi=0} as depending on the optimisation method it may be get stuck.
#' 
#' A default value for the tail fraction \code{phiu=TRUE} is given. 
#' The \code{\link[evmix:fnormgpd]{lnormgpd}} also has the usual defaults for
#' the other parameters, but \code{\link[evmix:fnormgpd]{nlnormgpd}} and
#' \code{\link[evmix:fnormgpd]{nlunormgpd}} has no defaults.
#' 
#' If the hessian is of reduced rank then the variance covariance (from inverse hessian)
#' and standard error of parameters cannot be calculated, then by default 
#' \code{std.err=TRUE} and the function will stop. If you want the parameter estimates
#' even if the hessian is of reduced rank (e.g. in a simulation study) then
#' set \code{std.err=FALSE}. 
#' 
#' Invalid parameter ranges will give \code{0} for likelihood, \code{log(0)=-Inf} for
#' log-likelihood and \code{-log(0)=Inf} for negative log-likelihood. 
#' 
#' Due to symmetry, the lower tail can be described by GPD by negating the data/quantiles.
#' 
#' Infinite and missing sample values are dropped.
#' 
#' Error checking of the inputs is carried out and will either stop or give warning message
#' as appropriate.
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
#' @section Acknowledgments: These functions are deliberately similar
#'   in syntax and functionality to the commonly used functions in the
#'   \code{\link[ismev:gpd.fit]{ismev}} and \code{\link[evd:fpot]{evd}} packages
#'   for which their author's contributions are gratefully acknowledged.
#'   
#'   Anna MacDonald and Xin Zhao laid some of the groundwork with programs they
#'   wrote for MATLAB.
#'   
#'   Clement Lee and Emma Eastoe suggested providing inbuilt
#'   profile likelihood estimation for threshold and fixed threshold approach.
#' 
#' @seealso \code{\link[stats:Normal]{dnorm}},
#'  \code{\link[evmix:fgpd]{fgpd}} and \code{\link[evmix:gpd]{gpd}}
#'
#' @aliases fnormgpd lnormgpd nlnormgpd proflunormgpd nlunormgpd
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
#' # Bulk model based tail fraction
#' fit = fnormgpd(x)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-4, 4))
#' lines(xx, y)
#' with(fit, lines(xx, dnormgpd(xx, nmean, nsd, u, sigmau, xi), col="red"))
#' abline(v = fit$u, col = "red")
#'   
#' # Parameterised tail fraction
#' fit2 = fnormgpd(x, phiu = FALSE)
#' with(fit2, lines(xx, dnormgpd(xx, nmean, nsd, u, sigmau, xi, phiu), col="blue"))
#' abline(v = fit2$u, col = "blue")
#' legend("topleft", c("True Density","Bulk Tail Fraction","Parameterised Tail Fraction"),
#'   col=c("black", "red", "blue"), lty = 1)
#'   
#' # Profile likelihood for initial value of threshold and fixed threshold approach
#' fitu = fnormgpd(x, useq = seq(0, 3, length = 20))
#' fitfix = fnormgpd(x, useq = seq(0, 3, length = 20), fixedu = TRUE)
#' 
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-4, 4))
#' lines(xx, y)
#' with(fit, lines(xx, dnormgpd(xx, nmean, nsd, u, sigmau, xi), col="red"))
#' abline(v = fit$u, col = "red")
#' with(fitu, lines(xx, dnormgpd(xx, nmean, nsd, u, sigmau, xi), col="purple"))
#' abline(v = fitu$u, col = "purple")
#' with(fitfix, lines(xx, dnormgpd(xx, nmean, nsd, u, sigmau, xi), col="darkgreen"))
#' abline(v = fitfix$u, col = "darkgreen")
#' legend("topleft", c("True Density","Default initial value (90% quantile)",
#'  "Prof. lik. for initial value", "Prof. lik. for fixed threshold"),
#'  col=c("black", "red", "purple", "darkgreen"), lty = 1)
#' }
#'   

# maximum likelihood fitting for normal bulk with GPD for upper tail
fnormgpd <- function(x, phiu = TRUE, useq = NULL, fixedu = FALSE, pvector = NULL,
  std.err = TRUE, method = "BFGS", control = list(maxit = 10000), finitelik = TRUE, ...) {

  call <- match.call()
    
  np = 5 # maximum number of parameters

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
      pvector[4] = initfgpd$sigmau
      pvector[5] = initfgpd$xi
    }
    
  } else { # profile or fixed
    
    check.nparam(pvector, nparam = np - 1, allownull = TRUE)

    # profile likelihood for threshold or scalar given
    if (length(useq) != 1) {
      
      # remove thresholds with less than 5 excesses
      useq = useq[sapply(useq, FUN = function(u, x) sum(x > u) > 5, x = x)]
      check.param(useq, allowvec = TRUE)
      
      nllhu = sapply(useq, proflunormgpd, pvector = pvector, x = x, phiu = phiu,
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
        pvector[3] = initfgpd$sigmau
        pvector[4] = initfgpd$xi
      }
    } else { # threshold free parameter
      if (is.null(pvector)) {
        pvector[1] = mean(x, trim = 0.2)
        pvector[2] = sd(x)
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
    nllh = nlunormgpd(pvector, u, x, phiu)
    if (is.infinite(nllh)) {
      pvector[4] = 0.1
      nllh = nlunormgpd(pvector, u, x, phiu)
    }
    if (is.infinite(nllh)) stop("initial parameter values are invalid")
  
    fit = optim(par = as.vector(pvector), fn = nlunormgpd, u = u, x = x, phiu = phiu,
      finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)    
    
    nmean = fit$par[1]
    nsd = fit$par[2]
    sigmau = fit$par[3]
    xi = fit$par[4]
    
  } else { # complete (non-separable) likelihood
    
    nllh = nlnormgpd(pvector, x, phiu)
    if (is.infinite(nllh)) {
      pvector[5] = 0.1
      nllh = nlnormgpd(pvector, x, phiu)
    }
    if (is.infinite(nllh)) stop("initial parameter values are invalid")
  
    fit = optim(par = as.vector(pvector), fn = nlnormgpd, x = x, phiu = phiu,
      finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)    
    
    nmean = fit$par[1]
    nsd = fit$par[2]
    u = fit$par[3]
    sigmau = fit$par[4]
    xi = fit$par[5]
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
#' @aliases fnormgpd lnormgpd nlnormgpd proflunormgpd nlunormgpd
#' @rdname  fnormgpd

# log-likelihood function for normal bulk with GPD for upper tail
lnormgpd <- function(x, nmean = 0, nsd = 1, u = qnorm(0.9, nmean, nsd),
  sigmau = nsd, xi = 0, phiu = TRUE, log = TRUE) {

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(nmean)
  check.param(nsd)
  check.param(u)
  check.param(sigmau)
  check.param(xi)
  check.phiu(phiu, allowfalse = TRUE)
  check.logic(log)

  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    x = x[is.finite(x)] # ignore missing and infinite cases
  }

  check.quant(x)
  n = length(x)

  check.inputn(c(length(nmean), length(nsd), length(u), length(sigmau), length(xi), length(phiu)),
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

  if ((nsd <= 0) | (sigmau <= 0) | (u <= min(x)) | (u >= max(x))) {
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
  
    syu = 1 + xi * (xu - u) / sigmau  
    yb = (xb - nmean) / nsd # used for normal
  
    if ((min(syu) <= 0) | (phiu <= 0) | (phiu >= 1) | (pu <= 0) | (pu >= 1)) {
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
#' @aliases fnormgpd lnormgpd nlnormgpd proflunormgpd nlunormgpd
#' @rdname  fnormgpd

# negative log-likelihood function for normal bulk with GPD for upper tail
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlnormgpd <- function(pvector, x, phiu = TRUE, finitelik = FALSE) {

  np = 5 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.phiu(phiu, allowfalse = TRUE)
  check.logic(finitelik)

  nmean = pvector[1]
  nsd = pvector[2]
  u = pvector[3]
  sigmau = pvector[4]
  xi = pvector[5]

  nllh = -lnormgpd(x, nmean, nsd, u, sigmau, xi, phiu) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}

#' @export
#' @aliases fnormgpd lnormgpd nlnormgpd proflunormgpd nlunormgpd
#' @rdname  fnormgpd

# profile negative log-likelihood function for given threshold for
# normal bulk with GPD for upper tail
# designed for sapply to loop over vector of thresholds (hence u is first input)
proflunormgpd <- function(u, pvector = NULL, x, phiu = TRUE, method = "BFGS",
  control = list(maxit = 10000), finitelik = TRUE, ...) {

  np = 5 # maximum number of parameters

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
    nllh = nlunormgpd(pvector, u, x, phiu)
    
    if (is.infinite(nllh)) pvector = NULL
  }

  if (is.null(pvector)) {
    pvector[1] = mean(x, trim = 0.2)
    pvector[2] = sd(x)
    initfgpd = fgpd(x, u, std.err = FALSE)
    pvector[3] = initfgpd$sigmau
    pvector[4] = initfgpd$xi
    nllh = nlunormgpd(pvector, u, x, phiu)
  }

  if (is.infinite(nllh)) {
    pvector[4] = 0.1
    nllh = nlunormgpd(pvector, u, x, phiu)
  }
  
  # if still invalid then output cleanly
  if (is.infinite(nllh)) {
    warning(paste("initial parameter values for threshold u =", u, "are invalid"))
    fit = list(par = rep(NA, np), value = Inf, counts = 0, convergence = NA, 
      message = "initial values invalid", hessian = rep(NA, np))
  } else {
    fit = optim(par = as.vector(pvector), fn = nlunormgpd, u = u, x = x, phiu = phiu,
      finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)
  }
    
  if (finitelik & is.infinite(fit$value)) {
    fit$value = sign(fit$value) * 1e6
  }

  fit$value
}

#' @export
#' @aliases fnormgpd lnormgpd nlnormgpd proflunormgpd nlunormgpd
#' @rdname  fnormgpd

# negative log-likelihood function for normal bulk with GPD for upper tail
# (wrapper for likelihood, designed for threshold to be fixed and other parameters optimised)
nlunormgpd <- function(pvector, u, x, phiu = TRUE, finitelik = FALSE) {

  np = 5 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np - 1)
  check.param(u)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.phiu(phiu, allowfalse = TRUE)
  check.logic(finitelik)
    
  nmean = pvector[1]
  nsd = pvector[2]
  sigmau = pvector[3]
  xi = pvector[4]

  nllh = -lnormgpd(x, nmean, nsd, u, sigmau, xi, phiu) 

  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}
