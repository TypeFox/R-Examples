#' @export
#' 
#' @title MLE Fitting of Generalised Pareto Distribution (GPD)
#'   
#' @description Maximum likelihood estimation for fitting the GPD with
#'   parameters scale \code{sigmau} and shape \code{xi} to the threshold
#'   exceedances, conditional on being above a threshold \code{u}. Unconditional
#'   likelihood fitting also provided when the probability \code{phiu} of being
#'   above the threshold \code{u} is given.
#'   
#' @param x          vector of sample data
#' @param u          scalar threshold
#' @param sigmau     scalar scale parameter (positive)
#' @param xi         scalar shape parameter
#' @param pvector    vector of initial values of GPD parameters (\code{sigmau}, \code{xi}) or \code{NULL}
#' @param phiu       probability of being above threshold \eqn{[0, 1]} or \code{NULL}, see Details
#' @param std.err    logical, should standard errors be calculated
#' @param method     optimisation method (see \code{\link[stats:optim]{optim}})
#' @param control    optimisation control list (see \code{\link[stats:optim]{optim}})
#' @param finitelik  logical, should log-likelihood return finite value for invalid parameters
#' @param log        logical, if \code{TRUE} then log-likelihood rather than likelihood is output
#' @param ...        optional inputs passed to \code{\link[stats:optim]{optim}}
#'   
#' @details The GPD is fitted to the exceedances of the threshold \code{u} using
#'   maximum likelihood estimation. The estimated parameters, 
#'   variance-covariance matrix and their standard errors are automatically 
#'   output.
#'   
#'   The log-likelihood and negative log-likelihood are also provided for wider 
#'   usage, e.g. constructing your own extreme value mixture model or profile
#'   likelihood functions. The 
#'   parameter vector \code{pvector} must be specified in the negative 
#'   log-likelihood \code{\link[evmix:fgpd]{nlgpd}}.
#'   
#'   Log-likelihood calculations are carried out in 
#'   \code{\link[evmix:fgpd]{lgpd}}, which takes parameters as inputs in the 
#'   same form as distribution functions. The negative log-likelihood is a 
#'   wrapper for \code{\link[evmix:fgpd]{lgpd}}, designed towards making it 
#'   useable for optimisation (e.g. parameters are given a vector as first 
#'   input).
#'   
#'   The default value for the tail fraction \code{phiu} in the fitting function
#'   \code{\link[evmix:fgpd]{fgpd}} is \code{NULL}, in which case the MLE is calculated 
#'   using the sample proportion of exceedances. In this case the standard error for \code{phiu} is 
#'   estimated and output as \code{se.phiu}, otherwise it is set to \code{NA}. Consistent with the 
#'   \code{\link[evd:fpot]{evd}} library the missing values (\code{NA} and 
#'   \code{NaN}) are assumed to be below the threshold in calculating the tail fraction.
#'   
#'   Otherwise, in the fitting function \code{\link[evmix:fgpd]{fgpd}} the tail 
#'   fraction \code{phiu} can be specified as any value over \eqn{(0, 1]}, i.e.
#'   excludes \eqn{\phi_u=0}, leading to the unconditional log-likelihood being
#'   used for estimation. In this case the standard error will be output as \code{NA}.
#'   
#'   In the log-likelihood functions \code{\link[evmix:fgpd]{lgpd}} and 
#'   \code{\link[evmix:fgpd]{nlgpd}} the tail fraction \code{phiu} cannot be
#'   \code{NULL} but can be over the range \eqn{[0, 1]}, i.e. which includes
#'   \eqn{\phi_u=0}.
#'   
#'   The value of \code{phiu} does not effect the GPD parameter estimates, only
#'   the value of the likelihood, as:
#'   
#'   \deqn{L(\sigma_u, \xi; u, \phi_u) = (\phi_u ^ {n_u}) L(\sigma_u, \xi; u,
#'   \phi_u=1)}
#'   
#'   where the GPD has scale \eqn{\sigma_u} and shape \eqn{\xi}, the threshold
#'   is \eqn{u} and \eqn{nu} is the number of exceedances. A non-unit value for
#'   \code{phiu} simply scales the likelihood and shifts the log-likelihood,
#'   thus the GPD parameter estimates are invariant to \code{phiu}.
#'   
#'   The default optimisation algorithm is "BFGS", which requires a finite
#'   negative log-likelihood function evaluation \code{finitelik=TRUE}. For
#'   invalid parameters, a zero likelihood is replaced with \code{exp(-1e6)}.
#'   The "BFGS" optimisation algorithms require finite values for likelihood, so
#'   any user input for \code{finitelik} will be overridden and set to
#'   \code{finitelik=TRUE} if either of these optimisation methods is chosen.
#'   
#'   It will display a warning for non-zero convergence result comes from 
#'   \code{\link[stats:optim]{optim}} function call.
#'   
#'   If the hessian is of reduced rank then the variance covariance (from
#'   inverse hessian) and standard error of parameters cannot be calculated,
#'   then by default \code{std.err=TRUE} and the function will stop. If you want
#'   the parameter estimates even if the hessian is of reduced rank (e.g. in a
#'   simulation study) then set \code{std.err=FALSE}.
#'   
#' @return \code{\link[evmix:fgpd]{lgpd}} gives (log-)likelihood and 
#'   \code{\link[evmix:fgpd]{nlgpd}} gives the negative log-likelihood. 
#'   \code{\link[evmix:fgpd]{fgpd}} returns a simple list with the following
#'   elements
#'   
#' \tabular{ll}{
#' \code{call}:     \tab \code{optim} call\cr
#' \code{x}:        \tab data vector \code{x}\cr
#' \code{init}:     \tab \code{pvector}\cr
#' \code{optim}:    \tab complete \code{optim} output\cr
#' \code{mle}:      \tab vector of MLE of parameters\cr
#' \code{cov}:      \tab variance-covariance matrix of MLE of parameters\cr
#' \code{se}:       \tab vector of standard errors of MLE of parameters\cr
#' \code{rate}:     \tab \code{phiu} to be consistent with \code{\link[evd:fpot]{evd}}\cr
#' \code{nllh}:     \tab minimum negative log-likelihood\cr
#' \code{n}:        \tab total sample size\cr
#' \code{u}:        \tab threshold\cr
#' \code{sigmau}:   \tab MLE of GPD scale\cr
#' \code{xi}:       \tab MLE of GPD shape\cr
#' \code{phiu}:     \tab MLE of tail fraction\cr
#' \code{se.phiu}:  \tab standard error of MLE of tail fraction (parameterised approach using sample proportion)\cr
#' }
#'   
#' The output list has some duplicate entries and repeats some of the inputs to both 
#' provide similar items to those from \code{\link[evd:fpot]{fpot}} and increase usability.
#'   
#' @note Unlike all the distribution functions for the GPD, the MLE fitting only
#'   permits single scalar values for each parameter, \code{phiu} and threshold
#'   \code{u}.
#'   
#'   When \code{pvector=NULL} then the initial values are calculated, type
#'   \code{fgpd} to see the default formulae used. The GPD fitting is not very
#'   sensitive to the initial values, so you will rarely have to  give
#'   alternatives. Avoid setting the starting value for the shape parameter to
#'   \code{xi=0} as depending on the optimisation method it may be get stuck.
#'   
#'   Default values for the threshold \code{u=0} and tail fraction
#'   \code{phiu=NULL} are given in the fitting \code{\link[evmix:fgpd]{fpgd}},
#'   in which case the MLE assumes that excesses over the threshold are given,
#'   rather than exceedances.
#'   
#'   The usual default of \code{phiu=1} is given in the likelihood functions
#'   \code{\link[evmix:fgpd]{lpgd}} and \code{\link[evmix:fgpd]{nlpgd}}.
#'   
#'   The \code{\link[evmix:fgpd]{lgpd}} also has the usual defaults for the
#'   other parameters, but \code{\link[evmix:fgpd]{nlgpd}} has no defaults.
#'   
#'   Infinite sample values are dropped in fitting function
#'   \code{\link[evmix:fgpd]{fpgd}}, but missing values are used to estimate
#'   \code{phiu} as described above. But in likelihood functions
#'   \code{\link[evmix:fgpd]{lpgd}} and \code{\link[evmix:fgpd]{nlpgd}} both
#'   infinite and missing values are ignored.
#'   
#'   Error checking of the inputs is carried out and will either stop or give
#'   warning message as appropriate.
#'   
#' @references
#' 
#' \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'   
#' @section Acknowledgments: Based on the \code{\link[ismev:gpd.fit]{gpd.fit}} and
#' \code{\link[evd:fpot]{fpot}} functions in the 
#' \code{\link[ismev:gpd.fit]{ismev}} and
#' \code{\link[evd:fpot]{evd}} packages for which their author's contributions are gratefully acknowledged.
#' They are designed to have similar syntax and functionality to simplify the transition for users of these packages.
#'   
#' @seealso \code{\link[evd:gpd]{dgpd}}, \code{\link[evd:fpot]{fpot}} and
#'   \code{\link[MASS:fitdistr]{fitdistr}}
#'   
#' @aliases fgpd lgpd nlgpd
#' @family  gpd fgpd
#'   
#' @examples
#' set.seed(1)
#' par(mfrow = c(2, 1))
#' 
#' # GPD is conditional model for threshold exceedances
#' # so tail fraction phiu not relevant when only have exceedances
#' x = rgpd(1000, u = 10, sigmau = 5, xi = 0.2)
#' xx = seq(0, 100, 0.1)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(0, 100))
#' lines(xx, dgpd(xx, u = 10, sigmau = 5, xi = 0.2))
#' fit = fgpd(x, u = 10)
#' lines(xx, dgpd(xx, u = fit$u, sigmau = fit$sigmau, xi = fit$xi), col="red")
#' 
#' # but tail fraction phiu is needed for conditional modelling of population tail
#' x = rnorm(10000)
#' xx = seq(-4, 4, 0.01)
#' hist(x, breaks = 200, freq = FALSE, xlim = c(0, 4))
#' lines(xx, dnorm(xx), lwd = 2)
#' fit = fgpd(x, u = 1)
#' lines(xx, dgpd(xx, u = fit$u, sigmau = fit$sigmau, xi = fit$xi, phiu = fit$phiu),
#'   col = "red", lwd = 2)
#' legend("topright", c("True Density","Fitted Density"), col=c("black", "red"), lty = 1)
#' 

# maximum likelihood fitting for GPD
fgpd <- function(x, u = 0, phiu = NULL, pvector = NULL, std.err = TRUE,
  method = "BFGS", control = list(maxit = 10000), finitelik = TRUE, ...) {

  call <- match.call()
    
  np = 2 # maximum number of parameters

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(u)
  check.prob(phiu, allownull = TRUE) # don't use check.phiu as TRUE only valid for mixture models
  check.nparam(pvector, nparam = np, allownull = TRUE)
  check.logic(std.err)
  check.optim(method)
  check.control(control)
  check.logic(finitelik)

  if (any(is.infinite(x))) warning("infinite cases have been removed")
  
  x = x[!is.infinite(x)] # ignore infinite cases only (all mixture models also ignore missing)
  
  if (any(is.na(x)))
    warning("missing values are treated as below threshold when estimating tail fraction")

  check.quant(x, allowna = TRUE)
  n = length(x)

  if ((method == "L-BFGS-B") | (method == "BFGS")) finitelik = TRUE

  # assume NA or NaN are below threshold consistent with evd library
  # hence use which() to ignore these

  xu = x[which(x > u)]
  nu = length(xu)
  
  if (nu < 1)
    stop("no elements of x are above threshold")
    
  # set default values if pvector is NULL
  if (is.null(pvector)) {
    yu = xu - u
    pvector[1] = sqrt(6 * var(yu)) / pi
    pvector[2] = 0.1
  }
  
  if (is.null(phiu)) {
    phiu = nu / n
    se.phiu = sqrt(phiu * (1 - phiu) / n)
  } else {
    if (phiu == 0) stop("tail probability must be in (0, 1]")
    se.phiu = NA
  }
  
  # check initial parameter vector and try different shape
  nllh = nlgpd(pvector, xu, u, phiu)
  if (is.infinite(nllh)) {
    pvector[2] = 0.1
    nllh = nlgpd(pvector, xu, u, phiu)
  }

  if (is.infinite(nllh)) stop("initial parameter values are invalid")

  fit = optim(par = as.vector(pvector), fn = nlgpd, x = xu, u = u, phiu = phiu,
    finitelik = finitelik, method = method, hessian = TRUE, ...)
  
  conv = TRUE
  if ((fit$convergence != 0) | any(fit$par == pvector) | (abs(fit$value) >= 1e6)) {
    conv = FALSE
    warning("check convergence")
  }
  
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
    n = n, u = u, sigmau = fit$par[1], xi = fit$par[2], phiu = phiu, se.phiu = se.phiu)
}

#' @export
#' @aliases fgpd lgpd nlgpd
#' @rdname  fgpd

# log-likelihood function for GPD
lgpd <- function(x, u = 0, sigmau = 1, xi = 0, phiu = 1, log = TRUE) {

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(u)
  check.param(sigmau)
  check.param(xi)
  check.prob(phiu) # don't use check.phiu as TRUE only valid for mixture models
  check.logic(log)

  check.inputn(c(length(u), length(sigmau), length(xi), length(phiu)), allowscalar = TRUE)

  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    x[!is.finite(x)] = NA # ignore missing and infinite cases
  }
  
  # assume NA or NaN are below threshold consistent with evd library
  # hence use which() to ignore these
  
  xu = x[which(x > u)]
  nu = length(xu)
  
  yu = (xu - u) / sigmau # used when shape is zero
  syu = 1 + xi * yu      # used when shape non-zero  
  
  if ((min(syu) <= 0) | (sigmau <= 0) | (phiu <= 0)  | (phiu > 1)) { # phiu = 1 is conditional likelihood
    l = -Inf
  } else {
    if (abs(xi) < 1e-6) {
      l = - nu * log(sigmau) - sum(yu) + nu * log(phiu)
    } else {
      l = - nu * log(sigmau) - (1/xi + 1) * sum(log(syu)) + nu * log(phiu)
    }
  }
  
  if (!log) l = exp(l)
  
  l
}

#' @export
#' @aliases fgpd lgpd nlgpd
#' @rdname  fgpd

# negative log-likelihood function for GPD
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlgpd <- function(pvector, x, u = 0, phiu = 1, finitelik = FALSE) {

  np = 2 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np) # must be provided
  check.param(u)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.prob(phiu) # don't use check.phiu as TRUE only valid for mixture models
  check.logic(finitelik)
  
  sigmau = pvector[1]
  xi = pvector[2]

  nllh = -lgpd(x, u, sigmau, xi, phiu) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }
  
  nllh
}
