#' @export
#' 
#' @title MLE Fitting of Hybrid Pareto Extreme Value Mixture Model
#'
#' @description Maximum likelihood estimation for fitting the hybrid Pareto extreme
#' value mixture model
#'
#' @param pvector vector of initial values of parameters
#'                (\code{nmean}, \code{nsd}, \code{xi}) or \code{NULL}
#' @inheritParams fnormgpd
#' @inheritParams dnormgpd
#' @inheritParams fgpd
#' 
#' @details The hybrid Pareto model is fitted to the entire dataset using maximum likelihood
#' estimation. The estimated parameters, variance-covariance matrix and their standard errors
#' are automatically output.
#' 
#' The log-likelihood and negative log-likelihood are also provided for wider
#' usage, e.g. constructing profile likelihood functions. The parameter vector
#' \code{pvector} must be specified in the negative log-likelihood
#' \code{\link[evmix:fhpd]{nlhpd}}.
#' 
#' Log-likelihood calculations are carried out in
#' \code{\link[evmix:fhpd]{lhpd}}, which takes parameters as inputs in
#' the same form as distribution functions. The negative log-likelihood is a
#' wrapper for \code{\link[evmix:fhpd]{lhpd}}, designed towards making
#' it useable for optimisation (e.g. parameters are given a vector as first
#' input).
#' 
#' Missing values (\code{NA} and \code{NaN}) are assumed to be invalid data so are ignored,
#' which is inconsistent with the \code{\link[evd:fpot]{evd}} library which assumes the 
#' missing values are below the threshold.
#'
#' The function \code{\link[evmix:fhpd]{lhpd}} carries out the calculations
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
#' \code{\link[stats:optim]{optim}} function call.
#' 
#' If the hessian is of reduced rank then the variance covariance (from inverse hessian)
#' and standard error of parameters cannot be calculated, then by default 
#' \code{std.err=TRUE} and the function will stop. If you want the parameter estimates
#' even if the hessian is of reduced rank (e.g. in a simulation study) then
#' set \code{std.err=FALSE}. 
#' 
#' @return \code{\link[evmix:fhpd]{lhpd}} gives (log-)likelihood and 
#' \code{\link[evmix:fhpd]{nlhpd}} gives the negative log-likelihood. 
#' \code{\link[evmix:fhpd]{fhpd}} returns a simple list with the following elements
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
#' \code{nmean}:  \tab MLE of normal mean\cr
#' \code{nsd}:    \tab MLE of normal standard deviation\cr
#' \code{u}:      \tab threshold\cr
#' \code{sigmau}: \tab MLE of GPD scale\cr
#' \code{xi}:     \tab MLE of GPD shape\cr
#' }
#' 
#' The output list has some duplicate entries and repeats some of the inputs to both 
#' provide similar items to those from \code{\link[evd:fpot]{fpot}} and to make it 
#' as useable as possible.
#'  
#' @note Unlike most of the distribution functions for the extreme value mixture models,
#' the MLE fitting only permits single scalar values for each parameter. Only the data is a vector.
#' 
#' When \code{pvector=NULL} then the initial values are calculated, type 
#' \code{fhpd} to see the default formulae used. The mixture model fitting can be
#' ***extremely*** sensitive to the initial values, so you if you get a poor fit then
#' try some alternatives. Avoid setting the starting value for the shape parameter to
#' \code{xi=0} as depending on the optimisation method it may be get stuck.
#' 
#' A default value for the tail fraction \code{phiu=TRUE} is given. 
#' The \code{\link[evmix:fhpd]{lhpd}} also has the usual defaults for
#' the other parameters, but \code{\link[evmix:fhpd]{nlhpd}} has no defaults.
#' 
#' Invalid parameter ranges will give \code{0} for likelihood, \code{log(0)=-Inf} for
#' log-likelihood and \code{-log(0)=Inf} for negative log-likelihood. 
#' 
#' Infinite and missing sample values are dropped.
#' 
#' Error checking of the inputs is carried out and will either stop or give warning message
#' as appropriate.
#' 
#' @references
#' \url{http://en.wikipedia.org/wiki/Normal_distribution}
#' 
#' \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
#' 
#' Scarrott, C.J. and MacDonald, A. (2012). A review of extreme value
#' threshold estimation and uncertainty quantification. REVSTAT - Statistical
#' Journal 10(1), 33-59. Available from \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}
#' 
#' Carreau, J. and Y. Bengio (2008). A hybrid Pareto model for asymmetric fat-tailed data:
#' the univariate case. Extremes 12 (1), 53-76.
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @seealso \code{\link[evmix:fgpd]{fgpd}} and \code{\link[evmix:gpd]{gpd}}
#' 
#' The \code{\link[condmixt:condmixt-package]{condmixt}} package written by one of the
#' original authors of the hybrid Pareto model (Carreau and Bengio, 2008) also has 
#' similar functions for the likelihood of the hybrid Pareto 
#' \code{\link[condmixt:hpareto.negloglike]{hpareto.negloglike}} and
#' fitting \code{\link[condmixt:hpareto.negloglike]{hpareto.fit}}.
#' 
#' @aliases fhpd lhpd nlhpd
#' @family  hpd hpdcon normgpd normgpdcon gng gngcon
#'          fhpd fhpdcon fnormgpd fnormgpdcon fgng fgngcon
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
#' # Hybrid Pareto provides reasonable fit for some asymmetric heavy upper tailed distributions
#' # but not for cases such as the normal distribution
#' fit = fhpd(x, std.err = FALSE)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-4, 4))
#' lines(xx, y)
#' with(fit, lines(xx, dhpd(xx, nmean, nsd, xi), col="red"))
#' abline(v = fit$u)
#' 
#' # Notice that if tail fraction is included a better fit is obtained
#' fit2 = fnormgpdcon(x, std.err = FALSE)
#' with(fit2, lines(xx, dnormgpdcon(xx, nmean, nsd, u, xi), col="blue"))
#' abline(v = fit2$u)
#' legend("topright", c("Standard Normal", "Hybrid Pareto", "Normal+GPD Continuous"),
#'   col=c("black", "red", "blue"), lty = 1)
#' }
#'  

# maximum likelihood fitting for hybrid Pareto
fhpd <- function(x, pvector = NULL, std.err = TRUE, method = "BFGS",
  control = list(maxit = 10000), finitelik = TRUE, ...) {
  
  call <- match.call()
  
  np = 3 # maximum number of parameters

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

  check.quant(x)
  n = length(x)

  if ((method == "L-BFGS-B") | (method == "BFGS")) finitelik = TRUE
  
  if (is.null(pvector)) {
    pvector[1] = mean(x)
    pvector[2] = sd(x)
    initfgpd = fgpd(x, as.vector(quantile(x, 0.9)), std.err = FALSE)
    pvector[3] = initfgpd$xi
  }
  
  nllh = nlhpd(pvector, x)
  if (is.infinite(nllh)) {
    pvector[3] = 0.1
    nllh = nlhpd(pvector, x)  
  }
  if (is.infinite(nllh)) stop("initial parameter values are invalid")

  fit = optim(par = as.vector(pvector), fn = nlhpd, x = x, finitelik = finitelik,
              method = method, control = control, hessian = TRUE, ...)
  
  conv = TRUE
  if ((fit$convergence != 0) | any(fit$par == pvector) | (abs(fit$value) >= 1e6)) {
    conv = FALSE
    warning("check convergence")
  }
  
  nmean = fit$par[1]
  nsd = fit$par[2]
  xi = fit$par[3]
  
  z = (1 + xi)^2/(2*pi)
  wz = lambert_W0(z)
  
  u = nmean + nsd * sqrt(wz) * sign(1 + xi)
  sigmau = nsd * abs(1 + xi) / sqrt(wz)

  r = 1 + pnorm(u, nmean, nsd)
  phiu = 1/r # same normalisation for GPD and normal

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
    conv = conv, cov = invhess, mle = fit$par, se = se, rate = phiu, nllh = fit$value,
    n = n, nmean = nmean, nsd = nsd, u = u, sigmau = sigmau, xi = xi, phiu = phiu)
}

#' @export
#' @aliases fhpd lhpd nlhpd
#' @rdname  fhpd

# log-likelihood function for hybrid Pareto
# will not stop evaluation unless it has to
lhpd <- function(x, nmean = 0, nsd = 1, xi = 0, log = TRUE) {
  
  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(nmean)
  check.param(nsd)
  check.param(xi)
  check.logic(log)

  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    x = x[is.finite(x)] # ignore missing and infinite cases
  }
  
  check.quant(x)
  n = length(x)
  
  check.inputn(c(length(nmean), length(nsd), length(xi)), allowscalar = TRUE)
  
  z = (1 + xi)^2/(2*pi)
  wz = lambert_W0(z)
  
  u = nmean + nsd * sqrt(wz) * sign(1 + xi)

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
    du = sqrt(wz)
    sigmau = nsd * abs(1 + xi) / du
    
    syu = 1 + xi * (xu - u) / sigmau  
    
    yb = (xb - nmean) / nsd # used for normal
    
    r = 1 + pnorm(u, nmean, nsd)

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
#' @aliases fhpd lhpd nlhpd
#' @rdname  fhpd

# negative log-likelihood function for hybrid Pareto extreme value mixture model
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlhpd <- function(pvector, x, finitelik = FALSE) {
  
  np = 3 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.logic(finitelik)
  
  nmean = pvector[1]
  nsd = pvector[2]
  xi = pvector[3]
    
  nllh = -lhpd(x, nmean, nsd, xi) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }
  
  nllh
}
