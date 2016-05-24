#' @export
#' 
#' @title MLE Fitting of Dynamically Weighted Mixture Model
#'
#' @description Maximum likelihood estimation for fitting the dynamically weighted mixture model
#'
#' @param pvector vector of initial values of parameters
#'                (\code{wshape}, \code{wscale}, \code{cmu}, \code{ctau}, \code{sigmau}, \code{xi}) or \code{NULL}
#' @inheritParams fnormgpd
#' @inheritParams dwm
#' @inheritParams fgpd
#' 
#' @details The dynamically weighted mixture model is fitted to the entire dataset using maximum 
#' likelihood estimation. The estimated parameters, variance-covariance matrix and their standard
#' errors are automatically output.
#' 
#' The log-likelihood and negative log-likelihood are also provided for wider
#' usage, e.g. constructing profile likelihood functions. The parameter vector
#' \code{pvector} must be specified in the negative log-likelihood \code{\link[evmix:fdwm]{nldwm}}.
#' 
#' Log-likelihood calculations are carried out in
#' \code{\link[evmix:fdwm]{ldwm}}, which takes parameters as inputs in
#' the same form as distribution functions. The negative log-likelihood is a
#' wrapper for \code{\link[evmix:fdwm]{ldwm}}, designed towards making
#' it useable for optimisation (e.g. parameters are given a vector as first
#' input).
#' 
#' Non-negative data are ignored.
#' 
#' Missing values (\code{NA} and \code{NaN}) are assumed to be invalid data so are ignored,
#' which is inconsistent with the \code{\link[evd:fpot]{evd}} library which assumes the 
#' missing values are below the threshold.
#' 
#' The default optimisation algorithm is "BFGS", which requires a finite negative 
#' log-likelihood function evaluation \code{finitelik=TRUE}. For invalid 
#' parameters, a zero likelihood is replaced with \code{exp(-1e6)}. The "BFGS" 
#' optimisation algorithms require finite values for likelihood, so any user 
#' input for \code{finitelik} will be overridden and set to \code{finitelik=TRUE} 
#' if either of these optimisation methods is chosen.
#' 
#' It will display a warning for non-zero convergence result comes from 
#' \code{\link[stats:optim]{optim}} function call.
#' 
#' If the hessian is of reduced rank then the variance covariance (from inverse hessian)
#' and standard error of parameters cannot be calculated, then by default 
#' \code{std.err=TRUE} and the function will stop. If you want the parameter estimates
#' even if the hessian is of reduced rank (e.g. in a simulation study) then
#' set \code{std.err=FALSE}. 
#' 
#' @return \code{\link[evmix:fdwm]{ldwm}} gives (log-)likelihood and 
#' \code{\link[evmix:fdwm]{nldwm}} gives the negative log-likelihood. 
#' \code{\link[evmix:fdwm]{fdwm}} returns a simple list with the following elements
#'
#' \tabular{ll}{
#' \code{call}:   \tab \code{optim} call\cr
#' \code{x}:      \tab data vector \code{x}\cr
#' \code{init}:   \tab \code{pvector}\cr
#' \code{optim}:  \tab complete \code{optim} output\cr
#' \code{mle}:    \tab vector of MLE of parameters\cr
#' \code{cov}:    \tab variance-covariance matrix of MLE of parameters\cr
#' \code{se}:     \tab vector of standard errors of MLE of parameters\cr
#' \code{rate}:   \tab \code{phiu} to be consistent with \code{\link[evd:fpot]{evd}}\cr
#' \code{nllh}:   \tab minimum negative log-likelihood\cr
#' \code{n}:      \tab total sample size\cr
#' \code{wshape}: \tab MLE of Weibull shape\cr
#' \code{wscale}: \tab MLE of Weibull scale\cr
#' \code{mu}:     \tab MLE of Cauchy location\cr
#' \code{tau}:    \tab MLE of Cauchy scale\cr
#' \code{sigmau}: \tab MLE of GPD scale\cr
#' \code{xi}:     \tab MLE of GPD shape\cr
#' }
#' 
#' The output list has some duplicate entries and repeats some of the inputs to both 
#' provide similar items to those from \code{\link[evd:fpot]{fpot}} and to make it 
#' as useable as possible.
#'  
#' @note Unlike most of the distribution functions for the extreme value mixture models,
#' the MLE fitting only permits single scalar values for each parameter and 
#' \code{phiu}. Only the data is a vector.
#' 
#' When \code{pvector=NULL} then the initial values are calculated, type 
#' \code{fdwm} to see the default formulae used. The mixture model fitting can be
#' ***extremely*** sensitive to the initial values, so you if you get a poor fit then
#' try some alternatives. Avoid setting the starting value for the shape parameter to
#' \code{xi=0} as depending on the optimisation method it may be get stuck.
#' 
#' Infinite and missing sample values are dropped.
#' 
#' Error checking of the inputs is carried out and will either stop or give warning message
#' as appropriate.
#' 
#' @references
#' \url{http://en.wikipedia.org/wiki/Weibull_distribution}
#' 
#' \url{http://en.wikipedia.org/wiki/Cauchy_distribution}
#' 
#' \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
#' 
#' Scarrott, C.J. and MacDonald, A. (2012). A review of extreme value
#' threshold estimation and uncertainty quantification. REVSTAT - Statistical
#' Journal 10(1), 33-59. Available from \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}
#' 
#' Frigessi, A., O. Haug, and H. Rue (2002). A dynamic mixture model for unsupervised tail
#' estimation without threshold selection. Extremes 5 (3), 219-235
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @section Acknowledgments: See Acknowledgments in
#'   \code{\link[evmix:fnormgpd]{fnormgpd}}, type \code{help fnormgpd}.
#'   
#' @seealso \code{\link[evmix:fgpd]{fgpd}} and \code{\link[evmix:gpd]{gpd}}
#' @aliases fdwm ldwm nldwm
#' @family  dwm fdwm
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' par(mfrow = c(1, 1))
#' 
#' x = rweibull(1000, shape = 2)
#' xx = seq(-0.1, 4, 0.01)
#' y = dweibull(xx, shape = 2)
#' 
#' fit = fdwm(x, std.err = FALSE)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-0.1, 4))
#' lines(xx, y)
#' with(fit, lines(xx, ddwm(xx, wshape, wscale, cmu, ctau, sigmau, xi), col="red"))
#' }
#' 

# maximum likelihood fitting for dynamically weighted mixture model
fdwm = function(x, pvector = NULL, std.err = TRUE, method = "BFGS",
  control = list(maxit = 10000), finitelik = TRUE, ...) {
  
  call <- match.call()
  
  np = 6 # maximum number of parameters
  
  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.nparam(pvector, nparam = np, allownull = TRUE)
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
    
  if (is.null(pvector)) {
    initfweibull = fitdistr(x, "weibull", lower = c(1e-8, 1e-8))
    pvector[1] = initfweibull$estimate[1]
    pvector[2] = initfweibull$estimate[2]    
    pvector[3] = quantile(x, 0.7)
    pvector[4] = sd(x)/10
    pvector[5] = sqrt(6*var(x))/pi
    pvector[6] = 0.1
  }
  
  nllh = nldwm(pvector, x)
  if (is.infinite(nllh)) stop("initial parameter values are invalid")

  fit = optim(par = as.vector(pvector), fn = nldwm, x = x, finitelik = finitelik,
    method = method, control = control, hessian = TRUE, ...)
  
  conv = TRUE
  if ((fit$convergence != 0) | any(fit$par == pvector) | (abs(fit$value) >= 1e6)) {
    conv = FALSE
    warning("check convergence")
  }
  
  wshape = fit$par[1]
  wscale = fit$par[2]
  cmu = fit$par[3]
  ctau = fit$par[4]
  sigmau = fit$par[5]
  xi = fit$par[6]
  
  if (conv & std.err) {
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
  
  list(call = call, x = as.vector(x), init = as.vector(pvector), optim = fit,
    conv = conv, cov = invhess, mle = fit$par, se = se, nllh = fit$value,
    n = n, wshape = wshape, wscale = wscale, cmu = cmu, ctau = ctau, sigmau = sigmau, xi = xi)
}

#' @export
#' @aliases fdwm ldwm nldwm
#' @rdname  fdwm

# log-likelihood function for dynamically weighted mixture model
# will not stop evaluation unless it has to
ldwm = function(x,  wshape = 1, wscale = 1, cmu = 1, ctau = 1,
  sigmau = sqrt(wscale^2 * gamma(1 + 2/wshape) - (wscale * gamma(1 + 1/wshape))^2),
  xi = 0, log = TRUE) {
  
  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(wshape)
  check.param(wscale)
  check.param(cmu)
  check.param(ctau)
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
  
  check.inputn(c(length(wshape), length(wscale), length(cmu), length(ctau), length(sigmau), length(xi)),
               allowscalar = TRUE)

  # assume NA or NaN are irrelevant as entire lower tail is now modelled
  # inconsistent with evd library definition
  # hence use which() to ignore these

  rx <- function(x, wshape, wscale, cmu, ctau, sigmau, xi) {
    (dgpd(x, 0, sigmau, xi) - dweibull(x, wshape, wscale))*atan((x - cmu)/ctau)
  }

  if ((wscale <= 0) | (wshape <= 0) | (cmu <= 0) | (ctau <= 0) | (sigmau <= 0) | (cmu >= max(x))) {
    l = -Inf
  } else {
        
    syu = 1 + xi * (x / sigmau) # zero threshold
    
    if (min(syu) <= 0) {
      l = -Inf
    } else { 
      px = pcauchy(x, cmu, ctau)

      r = try(integrate(rx, wshape, wscale, cmu = cmu, ctau = ctau, sigmau = sigmau, xi = xi,
        lower= 0, upper = Inf, subdivisions = 10000, rel.tol = 1.e-10, stop.on.error = FALSE)$value)
      
      if (inherits(r, "try-error")) {
        warning("numerical integration failed, ignore previous messages, optimisation will try again")
        l = -Inf
      } else {
        z = n*log((1 + r/pi))
        
        pweights = pcauchy(x, cmu, ctau)
        bulk = sum(log((1 - pweights)*dweibull(x, wshape, wscale) +
          pweights*dgpd(x, 0, sigmau, xi)))
        l = bulk-z
      }
    }
  }
  
  l
}

#' @export
#' @aliases fdwm ldwm nldwm
#' @rdname  fdwm

# negative log-likelihood function for dynamically weighted mixture model
# (wrapper for likelihood, inputs and checks designed for optimisation)
nldwm = function(pvector, x, finitelik = FALSE) {
  
  np = 6 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.logic(finitelik)
  
  wshape = pvector[1]
  wscale = pvector[2]
  cmu = pvector[3]
  ctau = pvector[4]
  sigmau = pvector[5]
  xi = pvector[6]
  
  nllh = -ldwm(x, wshape, wscale, cmu, ctau, sigmau, xi) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }
  
  nllh
}
