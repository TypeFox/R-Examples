#' @export
#' 
#' @title MLE Fitting of Mixture of Gammas Using EM Algorithm
#'
#' @description Maximum likelihood estimation for fitting the mixture of gammas distribution
#' using the EM algorithm.
#'
#' @param mgshape   mgamma shape (positive) as vector of length \code{M}
#' @param mgscale   mgamma scale (positive) as vector of length \code{M}
#' @param mgweight  mgamma weights (positive) as vector of length \code{M}
#' @param M         number of gamma components in mixture
#' @param tau       matrix of posterior probability of being in each component
#'                  (\code{nxM} where \code{n} is \code{length(x)})
#' @inheritParams   fgpd
#' 
#' @details The weighted mixture of gammas distribution is fitted to the entire
#' dataset by maximum likelihood estimation using the EM algorithm. The estimated parameters,
#' variance-covariance matrix and their standard errors are automatically output.
#' 
#' The expectation step estimates the expected probability of being in each component
#' conditional on gamma component parameters. The maximisation step optimizes the
#' negative log-likelihood conditional on posterior probabilities of each observation
#' being in each component.
#' 
#' The optimisation of the likelihood for these mixture models can be very sensitive to
#' the initial parameter vector, as often there are numerous local modes. This is an
#' inherent feature of such models and the EM algorithm. The EM algorithm is guaranteed
#' to reach the maximum of the local mode. Multiple initial values should be considered
#' to find the global maximum. If the \code{pvector} is input as \code{NULL} then 
#' random component probabilities are simulated as the initial values, so multiple such runs
#' should be run to check the sensitivity to initial values. Alternatives to black-box
#' likelihood optimisers (e.g. simulated annealing), or moving to computational Bayesian
#' inference, are also worth considering.
#' 
#' The log-likelihood functions are provided for wider usage, e.g. constructing profile
#' likelihood functions. The parameter vector \code{pvector} must be specified in the
#' negative log-likelihood functions \code{\link[evmix:fmgamma]{nlmgamma}} and
#' \code{\link[evmix:fmgamma]{nlEMmgamma}}.
#' 
#' Log-likelihood calculations are carried out in \code{\link[evmix:fmgamma]{lmgamma}},
#' which takes parameters as inputs in the same form as the distribution functions. The
#' negative log-likelihood function \code{\link[evmix:fmgamma]{nlmgamma}} is a wrapper
#' for \code{\link[evmix:fmgamma]{lmgamma}} designed towards making it useable for optimisation,
#' i.e. \code{\link[evmix:fmgamma]{nlmgamma}} has complete parameter vector as first input.
#' Similarly, for the maximisation step negative log-likelihood
#' \code{\link[evmix:fmgamma]{nlEMmgamma}}, which also has the second input as the component
#' probability vector \code{mgweight}.
#' 
#' Missing values (\code{NA} and \code{NaN}) are assumed to be invalid data so are ignored.
#' 
#' The function \code{\link[evmix:fnormgpd]{lnormgpd}} carries out the calculations
#' for the log-likelihood directly, which can be exponentiated to give actual
#' likelihood using (\code{log=FALSE}).
#' 
#' The default optimisation algorithm in the "maximisation step" is "BFGS", which
#' requires a finite negative 
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
#' Suppose there are \eqn{M} gamma components with (scalar) shape and scale parameters and
#' weight for each component. Only \eqn{M-1} are to be provided in the initial parameter
#' vector, as the \eqn{M}th components weight is uniquely determined from the others.
#' 
#' For the fitting function \code{\link[evmix:fmgamma]{fmgamma}} and negative log-likelihood
#' functions the parameter vector \code{pvector} is a \code{3*M-1} length vector
#' containing all \eqn{M} gamma component shape parameters first, 
#' followed by the corresponding \eqn{M} gamma scale parameters,
#' then all the corresponding \eqn{M-1} probability weight parameters. The full parameter vector
#' is then \code{c(mgshape, mgscale, mgweight[1:(M-1)])}.
#' 
#' For the maximisation step negative log-likelihood functions the parameter vector
#' \code{pvector} is a \code{2*M} length vector containing all \eqn{M} gamma component
#' shape parameters first followed by the corresponding \eqn{M} gamma scale parameters. The
#' partial parameter vector is then \code{c(mgshape, mgscale)}.
#'
#' For identifiability purposes the mean of each gamma component must be in ascending in order. 
#' If the initial parameter vector does not satisfy this constraint then an error is given. 
#' 
#' Non-positive data are ignored as likelihood is infinite, except for \code{gshape=1}.
#' 
#' @return Log-likelihood is given by \code{\link[evmix:fmgamma]{lmgamma}} and it's
#'   wrapper for negative log-likelihood from \code{\link[evmix:fmgamma]{nlmgamma}}. 
#'   The conditional negative log-likelihood
#'   using the posterior probabilities is given by \code{\link[evmix:fmgamma]{nlEMmgamma}}.
#'   Fitting function \code{\link[evmix:fmgammagpd]{fmgammagpd}} using EM algorithm returns
#'   a simple list with the following elements
#'   
#' \tabular{ll}{
#'  \code{call}:      \tab \code{optim} call\cr
#'  \code{x}:         \tab data vector \code{x}\cr
#'  \code{init}:      \tab \code{pvector}\cr
#'  \code{optim}:     \tab complete \code{optim} output\cr
#'  \code{mle}:       \tab vector of MLE of parameters\cr
#'  \code{cov}:       \tab variance-covariance matrix of MLE of parameters\cr
#'  \code{se}:        \tab vector of standard errors of MLE of parameters\cr
#'  \code{nllh}:      \tab minimum negative log-likelihood\cr
#'  \code{n}:         \tab total sample size\cr
#'  \code{M}:         \tab number of gamma components\cr
#'  \code{mgshape}:   \tab MLE of gamma shapes\cr
#'  \code{mgscale}:   \tab MLE of gamma scales\cr
#'  \code{mgweight}:  \tab MLE of gamma weights\cr
#'  \code{EMresults}: \tab EM results giving complete negative log-likelihood, estimated parameters
#'                         and conditional "maximisation step" negative log-likelihood and convergence result\cr
#'  \code{posterior}: \tab posterior probabilites\cr
#' }
#' 
#' @note In the fitting and profile likelihood functions, when \code{pvector=NULL} then the default initial values
#' are obtained under the following scheme:
#' \itemize{
#'  \item number of sample from each component is simulated from symmetric multinomial distribution;
#'  \item sample data is then sorted and split into groups of this size (works well when components
#'        have modes which are well separated);
#'  \item for data within each component approximate MLE's for the
#'        gamma shape and scale parameters are estimated.
#' }
#' The \code{\link[evmix:fmgamma]{lmgamma}}, \code{\link[evmix:fmgamma]{nlmgamma}} and
#' \code{\link[evmix:fmgamma]{nlEMmgamma}} have no defaults.
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
#' Infinite and missing sample values are dropped.
#' 
#' Error checking of the inputs is carried out and will either stop or give warning message
#' as appropriate.
#' 
#' @references
#' \url{http://www.math.canterbury.ac.nz/~c.scarrott/evmix}
#' 
#' \url{http://en.wikipedia.org/wiki/Gamma_distribution}
#' 
#' \url{http://en.wikipedia.org/wiki/Mixture_model}
#' 
#' McLachlan, G.J. and Peel, D. (2000). Finite Mixture Models. Wiley.
#' 
#' @author Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @section Acknowledgments: Thanks to Daniela Laas, University of St Gallen, Switzerland for reporting various bugs in these functions.
#' 
#' @seealso \code{\link[stats:GammaDist]{dgamma}} and \code{\link[mixtools:gammamixEM]{gammamixEM}}
#'  in \code{mixtools} package
#' 
#' @aliases fmgamma lmgamma nlmgamma nlEMmgamma
#' @family  mgamma fmgamma
#'          gammagpd gammagpdcon fgammagpd fgammagpdcon normgpd fnormgpd
#'          mgammagpd mgammagpdcon fmgammagpd fmgammagpdcon 
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' par(mfrow = c(1, 1))
#' 
#' x = c(rgamma(1000, shape = 1, scale = 1), rgamma(3000, shape = 6, scale = 2))
#' xx = seq(-1, 40, 0.01)
#' y = (dgamma(xx, shape = 1, scale = 1) + 3 * dgamma(xx, shape = 6, scale = 2))/4
#' 
#' # Fit by EM algorithm
#' fit = fmgamma(x, M = 2)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-1, 40))
#' lines(xx, y)
#' with(fit, lines(xx, dmgamma(xx, mgshape, mgscale, mgweight), col="red"))
#' }
#'

# maximum likelihood fitting for mixture of gammas bulk with GPD for upper tail
fmgamma <- function(x, M, pvector = NULL, std.err = TRUE,
  method = "BFGS", control = list(maxit = 10000), finitelik = TRUE, ...) {

  call <- match.call()
    
  check.n(M)
  if (M == 1) stop("use fitdistr instead")
  np = 3*M # maximum number of parameters

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
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
  
  if (is.null(pvector)) {
    if (is.unsorted(x)) {
      x = sort(x)
    } else {
      if (x[1] > x[length(x)])
        x = rev(x)
    }
    
    # split dataset into components
    nM = rmultinom(1, n, rep(1/M, M))

    # try further simulations if any component is too small
    i = 1
    maxi = 10
    while (any(nM <= 2) & (i < maxi)) {
      nM = rmultinom(1, n, rep(1/M, M))
      i = i + 1
    }
     
    if (i == maxi) stop("sample is too small, user specified initial value (pvector) needed")
    
    # approximate MLE for each component
    approxmle = function(x) {
      s = log(mean(x)) - mean(log(x))
      k = (3 - s + sqrt((s - 3) ^ 2 + 24 * s))/12/s
      c(k, mean(x)/k)
    }
    
    mgparam = simplify2array(tapply(x, INDEX = rep(1:M, times = nM), FUN = approxmle, simplify=TRUE))

    pvector = c(mgparam[1,], mgparam[2,], nM[-M]/n)
  }

  if (length(pvector) != (np - 1)) 
    stop(paste("initial parameter vector must be of length", np - 1))

  check.nparam(pvector, nparam = np - 1, allownull = FALSE)

  if ((method == "L-BFGS-B") | (method == "BFGS")) finitelik = TRUE
  
  # EM algorithm convergence conditions and constraints
  maxit = 1000
  abstol = 1e-8
  
  # Initial values
  mgshape = pvector[1:M]
  mgscale = pvector[(M + 1):(2*M)]
  mgweight = pvector[(2*M + 1):(3*M - 1)]
  mgweight = c(mgweight, 1 - sum(mgweight)) # add constrained weight

  # check ordering of means, required for identifiability
  mui = mgshape * mgscale

  if (any(diff(mui) <= 0))
    stop("initial parameter vector does not satisfy constraint of strictly increasing component means")

  # Store iteration results
  nllh = rep(NA, maxit + 1)
  nllh[1] = 0
  
  EMresults = as.data.frame(matrix(NA, nrow = maxit, ncol = 3*M + 3))
  names(EMresults) = c("nllh", paste("mgshape", 1:M, sep = ""), paste("mgscale", 1:M, sep = ""),
    paste("mgweight", 1:M, sep = ""), "nllh.cond.weights", "conv.cond.weights")
  
  # approximate negative log-likelihood at initial step
  i = 1
  nllh[i + 1] = nlmgamma(pvector, x, M)
  
  if (is.infinite(nllh[i + 1])) stop("initial parameter values are invalid")

  EMresults[i,] = c(nllh[i + 1], mgshape, mgscale, mgweight, NA, NA)

  while ((abs(nllh[i + 1] - nllh[i]) > abstol) & (i < maxit)) {

    i = i + 1

    # Expectation step
    Mshape = matrix(mgshape, nrow = n, ncol = M, byrow = TRUE)
    Mscale = matrix(mgscale, nrow = n, ncol = M, byrow = TRUE)

    Mprob = t(mgweight * t(dgamma(x, Mshape, scale = Mscale)))
    
    # Estimated (posterior) probability of individual data being from each component,
    # conditional on each components' parameters
    tau = Mprob / rowSums(Mprob)
   
    # Updated mixture probabilities
    mgweight = colMeans(tau)
    
    # Maximisation step
    # MLE of each components' parameters,
    # conditional on probability of coming from each each component from E-step
  
    fit = try(optim(par = c(mgshape, mgscale), fn = nlEMmgamma,
      tau = tau, mgweight = mgweight, x = x, M = M,
      finitelik = finitelik, method = method, control = control, hessian = TRUE, ...))

    if (inherits(fit, "try-error")) stop("Maximisation step failed, try new initial values")

    conv = TRUE
    if ((fit$convergence != 0) | (abs(fit$value) >= 1e6)) {
      conv = FALSE
      warning("check convergence")
    }

    # Updated mixture component parameters
    mgshape = fit$par[1:M]
    mgscale = fit$par[(M + 1):(2*M)]

    # store results
    nllh[i + 1] = nlmgamma(c(mgshape, mgscale, mgweight[-M]), x, M)    
    EMresults[i,] = c(nllh[i], mgshape, mgscale, mgweight, fit$value, fit$convergence)
  }
  
  if (i == maxit) {
    warning("Maximum number of iteration reached, could be that M is too large or poor starting values")
  } else {
    print(paste(i, "iterations of EM algorithm"))
  }

  # Only output actual results, not blank rows
  EMresults = EMresults[1:i,]
  
  if ((fit$convergence != 0) | any(fit$par == pvector[1:(2*M)]) |
      any(mgweight[-M] == pvector[(2*M + 1):(3*M - 1)]) | (abs(fit$value) >= 1e6)) {
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
  
  list(call = call, x = as.vector(x), init = as.vector(pvector),
    optim = fit, conv = conv, cov = invhess, mle = fit$par, se = se,
    nllh = nllh[i + 1], n = n, M = M, mgshape = mgshape, mgscale = mgscale, mgweight = mgweight,
    EMresults = EMresults, posterior = tau)
}

#' @export
#' @aliases fmgamma lmgamma nlmgamma nlEMmgamma
#' @rdname  fmgamma

# log-likelihood function for mixture of gammas
lmgamma <- function(x, mgshape, mgscale, mgweight, log = TRUE) {
  
  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.logic(log)

  # user may try to input lists for mixture of gammas parameter vectors
  if (is.list(mgshape)) mgshape = unlist(mgshape)
  if (is.list(mgscale)) mgscale = unlist(mgscale)
  if (is.list(mgweight)) mgweight = unlist(mgweight)

  check.param(mgshape, allowvec = TRUE)
  check.param(mgscale, allowvec = TRUE)
  check.param(mgweight, allowvec = TRUE)

  # How many components in mixture
  M = check.inputn(c(length(mgshape), length(mgshape), length(mgweight)), allowscalar = TRUE)

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

  # gamma components must be ordered to ensure identifiability (means strictly ascending order)
  mui = mgshape * mgscale

  if (any(mgscale <= 0) | any(mgshape <= 0) | any(diff(mui) <= 0) |
      any(mgweight <= 0) | (sum(mgweight) > 1)) {
    l = -Inf
  } else {
    l = sum(dmgamma(x, mgshape, mgscale, mgweight, log = TRUE))
  }
  
  if (!log) l = exp(l)
  
  l
}

#' @export
#' @aliases fmgamma lmgamma nlmgamma nlEMmgamma
#' @rdname  fmgamma

# negative log-likelihood function for mixture of gammas
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlmgamma <- function(pvector, x, M, finitelik = FALSE) {

  check.n(M)
  if (M == 1) stop("use dgamma instead")
  np = 3*M # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np - 1)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.logic(finitelik)

  mgshape = pvector[1:M]
  mgscale = pvector[(M + 1):(2*M)]
  mgweight = pvector[(2*M + 1):(3*M - 1)]
  mgweight = c(mgweight, 1 - sum(mgweight)) # add constrained weight    

  nllh = -lmgamma(x, mgshape, mgscale, mgweight) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}

#' @export
#' @aliases fmgamma lmgamma nlmgamma nlEMmgamma
#' @rdname  fmgamma

# negative log-likelihood function for mixture of gammas,
# with component probabilities separated for EM algorithm
nlEMmgamma <- function(pvector, tau, mgweight, x, M, finitelik = FALSE) {

  check.n(M)
  if (M == 1) stop("use dgamma instead")
  np.noweight = 2*M # maximum number of non-weight parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np.noweight)
  check.nparam(mgweight, nparam = M)
  check.prob(mgweight)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
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

  if (!is.matrix(tau)) stop("Posterior probabilities (tau) must be nxM probability matrix")
  if (any(dim(tau) != c(n, M))) stop("Posterior probabilities (tau) must be nxM probability matrix")
  if (any((tau < 0) | (tau > 1))) stop("Posterior probabilities (tau) must be nxM probability matrix")
  
  mgshape = pvector[1:M]
  mgscale = pvector[(M + 1):(2*M)]
  
  # gamma components must be ordered to ensure identifiability (means strictly ascending order)
  mui = mgshape * mgscale

  if (any(mgscale <= 0) | any(mgshape <= 0) | any(diff(mui) <= 0)) {
    nllh = Inf
  } else {
    Mshape = matrix(mgshape, nrow = n, ncol = M, byrow = TRUE)
    Mscale = matrix(mgscale, nrow = n, ncol = M, byrow = TRUE)
  
    nllh = -sum(tau * log(t(mgweight * t(dgamma(x, Mshape, scale = Mscale)))))
  }
    
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}
