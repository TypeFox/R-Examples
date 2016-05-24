#' @export
#' 
#' @title MLE Fitting of Mixture of Gammas Bulk and GPD Tail Extreme Value Mixture Model
#' with Single Continuity Constraint using the EM algorithm.
#'
#' @description Maximum likelihood estimation for fitting the extreme value 
#' mixture model with mixture of gammas for bulk distribution upto the threshold and conditional
#' GPD above threshold with continuity at threshold. With options for profile likelihood estimation for threshold and
#' fixed threshold approach.
#'
#' @param mgshape   mgamma shape (positive) as vector of length \code{M}
#' @param mgscale   mgamma scale (positive) as vector of length \code{M}
#' @param mgweight  mgamma weights (positive) as vector of length \code{M}
#' @param M         number of gamma components in mixture
#' @param tau       matrix of posterior probability of being in each component
#'                  (\code{nxM} where \code{n} is \code{length(x)})
#' @inheritParams fnormgpd
#' @inheritParams fgpd
#' 
#' @details The extreme value mixture model with weighted mixture of gammas bulk and GPD tail with continuity at threshold is 
#' fitted to the entire dataset using maximum likelihood estimation using the EM algorithm. The estimated
#' parameters, variance-covariance matrix and their standard errors are automatically
#' output.
#' 
#' See help for \code{\link[evmix:fnormgpd]{fnormgpd}} for details, type \code{help fnormgpd}. 
#' Only the different features are outlined below for brevity.
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
#' negative log-likelihood functions \code{\link[evmix:fmgammagpdcon]{nlmgammagpdcon}} and
#' \code{\link[evmix:fmgammagpdcon]{nlEMmgammagpdcon}}.
#' 
#' Log-likelihood calculations are carried out in \code{\link[evmix:fmgammagpdcon]{lmgammagpdcon}},
#' which takes parameters as inputs in the same form as the distribution functions. The
#' negative log-likelihood function \code{\link[evmix:fmgammagpdcon]{nlmgammagpdcon}} is a wrapper
#' for \code{\link[evmix:fmgammagpdcon]{lmgammagpdcon}} designed towards making it useable for optimisation,
#' i.e. \code{\link[evmix:fmgammagpdcon]{nlmgammagpdcon}} has complete parameter vector as first input.
#' Though it is not directly used for optimisation here, as the EM algorithm due to mixture of
#' gammas for the bulk component of this model
#' 
#' The EM algorithm for the mixture of gammas utilises the
#' negative log-likelihood function \code{\link[evmix:fmgammagpdcon]{nlEMmgammagpdcon}}
#' which takes the posterior probabilities \eqn{tau} and component probabilities
#' \code{mgweight} as secondary inputs.
#' 
#' The profile likelihood for the threshold \code{\link[evmix:fmgammagpdcon]{proflumgammagpdcon}}
#' also implements the EM algorithm for the mixture of gammas, utilising the negative
#' log-likelihood function \code{\link[evmix:fmgammagpdcon]{nluEMmgammagpdcon}} which takes
#' the threshold, posterior probabilities \eqn{tau} and component probabilities
#' \code{mgweight} as secondary inputs. 
#' 
#' Missing values (\code{NA} and \code{NaN}) are assumed to be invalid data so are ignored.
#' 
#' Suppose there are \eqn{M} gamma components with (scalar) shape and scale parameters and
#' weight for each component. Only \eqn{M-1} are to be provided in the initial parameter
#' vector, as the \eqn{M}th components weight is uniquely determined from the others.
#' 
#' The initial parameter vector \code{pvector} always has the \eqn{M} gamma component
#' shape parameters followed by the corresponding \eqn{M} gamma scale parameters. However,
#' subsets of the other parameters are needed depending on which function is being used:
#' \itemize{
#'  \item {fmgammagpdcon} - \code{c(mgshape, mgscale, mgweight[1:(M-1)], u, xi)}
#'  \item {nlmgammagpdcon} - \code{c(mgshape, mgscale, mgweight[1:(M-1)], u, xi)}
#'  \item {nlumgammagpdcon} and {proflumgammagpdcon} - \code{c(mgshape, mgscale, mgweight[1:(M-1)], xi)}
#'  \item {nlEMmgammagpdcon} - \code{c(mgshape, mgscale, u, xi)}
#'  \item {nluEMmgammagpdcon} - \code{c(mgshape, mgscale, xi)}
#' }
#' Notice that when the component probability weights are included only the first \eqn{M-1} 
#' are specified, as the remaining one can be uniquely determined from these. Where some
#' parameters are left out, they are always taken as secondary inputs to the functions.
#' 
#' For identifiability purposes the mean of each gamma component must be in ascending in order. 
#' If the initial parameter vector does not satisfy this constraint then an error is given. 
#' 
#' Non-positive data are ignored as likelihood is infinite, except for \code{gshape=1}.
#' 
#' @return Log-likelihood is given by \code{\link[evmix:fmgammagpdcon]{lmgammagpdcon}} and it's
#'   wrappers for negative log-likelihood from \code{\link[evmix:fmgammagpdcon]{nlmgammagpdcon}}
#'   and \code{\link[evmix:fmgammagpdcon]{nlumgammagpdcon}}. The conditional negative log-likelihoods
#'   using the posterior probabilities are  \code{\link[evmix:fmgammagpdcon]{nlEMmgammagpdcon}}
#'   and \code{\link[evmix:fmgammagpdcon]{nluEMmgammagpdcon}}. Profile likelihood for single
#'   threshold given by \code{\link[evmix:fmgammagpdcon]{proflumgammagpdcon}} using EM algorithm. Fitting function
#'   \code{\link[evmix:fmgammagpdcon]{fmgammagpdcon}} using EM algorithm returns a simple list with the
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
#'  \code{M}:         \tab number of gamma components\cr
#'  \code{mgshape}:   \tab MLE of gamma shapes\cr
#'  \code{mgscale}:   \tab MLE of gamma scales\cr
#'  \code{mgweight}:  \tab MLE of gamma weights\cr
#'  \code{u}:         \tab threshold (fixed or MLE)\cr
#'  \code{sigmau}:    \tab MLE of GPD scale\cr
#'  \code{xi}:        \tab MLE of GPD shape\cr
#'  \code{phiu}:      \tab MLE of tail fraction (bulk model or parameterised approach)\cr
#'  \code{se.phiu}:   \tab standard error of MLE of tail fraction\cr
#'  \code{EMresults}: \tab EM results giving complete negative log-likelihood, estimated parameters
#'                         and conditional "maximisation step" negative log-likelihood and convergence result\cr
#'  \code{posterior}: \tab posterior probabilites\cr
#' }
#' 
#' @note In the fitting and profile likelihood functions, when \code{pvector=NULL} then the
#' default initial values are obtained under the following scheme:
#' \itemize{
#'  \item number of sample from each component is simulated from symmetric multinomial distribution;
#'  \item sample data is then sorted and split into groups of this size (works well when components
#'        have modes which are well separated);
#'  \item for data within each component approximate MLE's for the
#'        gamma shape and scale parameters are estimated;
#'  \item threshold is specified as sample 90\% quantile; and 
#'  \item MLE of GPD shape parameter above threshold. 
#' }
#' The other likelihood functions \code{\link[evmix:fmgammagpdcon]{lmgammagpdcon}},
#' \code{\link[evmix:fmgammagpdcon]{nlmgammagpdcon}}, \code{\link[evmix:fmgammagpdcon]{nlumgammagpdcon}} and
#' \code{\link[evmix:fmgammagpdcon]{nlEMmgammagpdcon}} and \code{\link[evmix:fmgammagpdcon]{nluEMmgammagpdcon}}
#' have no defaults.
#' 
#' @references
#' \url{http://www.math.canterbury.ac.nz/~c.scarrott/evmix}
#' 
#' \url{http://en.wikipedia.org/wiki/Gamma_distribution}
#' 
#'  \url{http://en.wikipedia.org/wiki/Mixture_model}
#' 
#' \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
#' 
#' McLachlan, G.J. and Peel, D. (2000). Finite Mixture Models. Wiley.
#' 
#' Scarrott, C.J. and MacDonald, A. (2012). A review of extreme value
#' threshold estimation and uncertainty quantification. REVSTAT - Statistical
#' Journal 10(1), 33-59. Available from \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}
#' 
#' Hu, Y. (2013). Extreme value mixture modelling: An R package and simulation study.
#' MSc (Hons) thesis, University of Canterbury, New Zealand.
#' \url{http://ir.canterbury.ac.nz/simple-search?query=extreme&submit=Go}
#' 
#' do Nascimento, F.F., Gamerman, D. and Lopes, H.F. (2011). A semiparametric
#' Bayesian approach to extreme value estimation. Statistical Computing, 22(2), 661-675.
#' 
#' @author Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @section Acknowledgments: Thanks to Daniela Laas, University of St Gallen, Switzerland for reporting various bugs in these functions.
#' 
#' See Acknowledgments in
#'   \code{\link[evmix:fnormgpd]{fnormgpd}}, type \code{help fnormgpd}.
#' 
#' @seealso \code{\link[stats:GammaDist]{dgamma}},
#'  \code{\link[evmix:fgpd]{fgpd}} and \code{\link[evmix:gpd]{gpd}}
#'  
#' @aliases fmgammagpdcon lmgammagpdcon nlmgammagpdcon nlumgammagpdcon proflumgammagpdcon nlEMmgammagpdcon nluEMmgammagpdcon
#' @family  mgamma fmgamma
#'          gammagpd gammagpdcon fgammagpd fgammagpdcon normgpd fnormgpd
#'          mgammagpd mgammagpdcon fmgammagpd fmgammagpdcon 
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' par(mfrow = c(2, 1))
#' 
#' n=1000
#' x = c(rgamma(n*0.25, shape = 1, scale = 1), rgamma(n*0.75, shape = 6, scale = 2))
#' xx = seq(-1, 40, 0.01)
#' y = (0.25*dgamma(xx, shape = 1, scale = 1) + 0.75 * dgamma(xx, shape = 6, scale = 2))
#' 
#' # Bulk model based tail fraction
#' # very sensitive to initial values, so best to provide sensible ones
#' fit.noinit = fmgammagpdcon(x, M = 2)
#' fit.withinit = fmgammagpdcon(x, M = 2, pvector = c(1, 6, 1, 2, 0.5, 15, 0.1))
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-1, 40))
#' lines(xx, y)
#' with(fit.noinit, lines(xx, dmgammagpdcon(xx, mgshape, mgscale, mgweight, u, xi), col="red"))
#' abline(v = fit.noinit$u, col = "red")
#' with(fit.withinit, lines(xx, dmgammagpdcon(xx, mgshape, mgscale, mgweight, u, xi), col="green"))
#' abline(v = fit.withinit$u, col = "green")
#'   
#' # Parameterised tail fraction
#' fit2 = fmgammagpdcon(x, M = 2, phiu = FALSE, pvector = c(1, 6, 1, 2, 0.5, 15, 0.1))
#' with(fit2, lines(xx, dmgammagpdcon(xx, mgshape, mgscale, mgweight, u, xi, phiu), col="blue"))
#' abline(v = fit2$u, col = "blue")
#' legend("topright", c("True Density","Default pvector", "Sensible pvector",
#'  "Parameterised Tail Fraction"), col=c("black", "red", "green", "blue"), lty = 1)
#'   
#' # Fixed threshold approach
#' fitfix = fmgammagpdcon(x, M = 2, useq = 15, fixedu = TRUE,
#'    pvector = c(1, 6, 1, 2, 0.5, 0.1))
#' 
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-1, 40))
#' lines(xx, y)
#' with(fit.withinit, lines(xx, dmgammagpdcon(xx, mgshape, mgscale, mgweight, u, xi), col="red"))
#' abline(v = fit.withinit$u, col = "red")
#' with(fitfix, lines(xx, dmgammagpdcon(xx,mgshape, mgscale, mgweight, u, xi), col="darkgreen"))
#' abline(v = fitfix$u, col = "darkgreen")
#' legend("topright", c("True Density", "Default initial value (90% quantile)",
#'  "Fixed threshold approach"), col=c("black", "red", "darkgreen"), lty = 1)
#' }
#'   

# maximum likelihood fitting for mixture of gammas bulk with GPD for upper tail with continuity at threshold
fmgammagpdcon <- function(x, M, phiu = TRUE, useq = NULL, fixedu = FALSE, pvector = NULL,
  std.err = TRUE, method = "BFGS", control = list(maxit = 10000), finitelik = TRUE, ...) {

  call <- match.call()
    
  check.n(M)
  if (M == 1) stop("use fgammagpdcon instead")
  np = 3*M + 2 # maximum number of parameters

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

  if (is.unsorted(x)) {
    x = sort(x)
  } else {
    if (x[1] > x[length(x)])
      x = rev(x)
  }

  if ((method == "L-BFGS-B") | (method == "BFGS")) finitelik = TRUE
  
  # useq must be specified if threshold is fixed
  if (fixedu & is.null(useq))
    stop("for fixed threshold approach, useq must be specified (as scalar or vector)")
  
  # approximate MLE for each gamma component
  approxmle = function(x) {
    s = log(mean(x)) - mean(log(x))
    k = (3 - s + sqrt((s - 3) ^ 2 + 24 * s))/12/s
    c(k, mean(x)/k)
  }

  # Check if profile likelihood or fixed threshold is being used
  # and determine initial values for parameters in each case
  if (is.null(useq)) { # not profile or fixed

    check.nparam(pvector, nparam = np - 1, allownull = TRUE)
    
    if (is.null(pvector)) {
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
      
      mgparam = simplify2array(tapply(x, INDEX = rep(1:M, times = nM), FUN = approxmle, simplify=TRUE))

      pvector = c(mgparam[1,], mgparam[2,], nM[-M]/n)
    
      # GPD defaults
      pvector[3*M] = as.vector(quantile(x, 0.9))
      initfgpd = fgpd(x, pvector[3*M], std.err = FALSE)
      pvector[3*M + 1] = initfgpd$xi
    }
    
  } else { # profile or fixed
    
    check.nparam(pvector, nparam = np - 2, allownull = TRUE)

    # profile likelihood for threshold or scalar given
    if (length(useq) != 1) {
      
      # remove thresholds with less than 5 excesses
      useq = useq[sapply(useq, FUN = function(u, x) sum(x > u) > 5, x = x)]
      check.posparam(useq, allowvec = TRUE)
      
      nllhu = sapply(useq, proflumgammagpdcon, pvector = pvector, x = x, M = M, phiu = phiu,
        method = method, control = control, finitelik = finitelik, ...)
      
      if (all(!is.finite(nllhu))) stop("thresholds are all invalid")
      u = useq[which.min(nllhu)]

    } else {
      u = useq
    }

    if (fixedu) { # threshold fixed
      if (is.null(pvector)) {
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
      
        mgparam = simplify2array(tapply(x, INDEX = rep(1:M, times = nM), FUN = approxmle, simplify=TRUE))

        pvector = c(mgparam[1,], mgparam[2,], nM[-M]/n)
    
        # GPD defaults
        initfgpd = fgpd(x, u, std.err = FALSE)
        pvector[3*M] = initfgpd$xi
      }
    } else { # threshold as initial value in usual MLE
      if (is.null(pvector)) {
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
      
        mgparam = simplify2array(tapply(x, INDEX = rep(1:M, times = nM), FUN = approxmle, simplify=TRUE))

        pvector = c(mgparam[1,], mgparam[2,], nM[-M]/n)
    
        # GPD defaults
        pvector[3*M] = u
        initfgpd = fgpd(x, u, std.err = FALSE)
        pvector[3*M + 1] = initfgpd$xi
      } else {
        pvector[3*M + 1] = pvector[3*M] # shift GPD shape to add in u
        pvector[3*M] = u
      }
    }
  }

  check.param(pvector, allowvec = TRUE)

  # EM algorithm convergence conditions and constraints
  maxit = 1000
  abstol = 1e-8

  # Note each component is treated as gamma+gpd, where gpd is same in each but scaling
  # is phiu * mgweight and each gamma is scaled by phib * mgweight

  if (fixedu) { # fixed threshold (separable) likelihood
    
    # Initial values
    mgshape = pvector[1:M]
    mgscale = pvector[(M + 1):(2*M)]
    mgweight = pvector[(2*M + 1):(3*M - 1)]
    mgweight = c(mgweight, 1 - sum(mgweight)) # add constrained weight
    xi = pvector[3*M]
    
    # check ordering of means, required for identifiability
    mui = mgshape * mgscale
    
    if (any(diff(mui) <= 0))
      stop("initial parameter vector does not satisfy constraint of strictly increasing component means")
    
    # Store iteration results
    nllh = rep(NA, maxit + 1)
    nllh[1] = 0
    
    EMresults = as.data.frame(matrix(NA, nrow = maxit, ncol = 3*M + 2 + 2))
    names(EMresults) = c("nllh", paste("mgshape", 1:M, sep = ""), paste("mgscale", 1:M, sep = ""),
      paste("mgweight", 1:M, sep = ""), "xi", "nllh.cond.weights", "conv.cond.weights")
    
    # approximate negative log-likelihood at initial step
    i = 1
    nllh[i + 1] = nlumgammagpdcon(pvector, u, x, M, phiu)  
    
    if (is.infinite(nllh[i + 1])) stop("initial parameter values are invalid")
    
    EMresults[i,] = c(nllh[i + 1], mgshape, mgscale, mgweight, xi, NA, NA)
    
    while ((abs(nllh[i + 1] - nllh[i]) > abstol) & (i < maxit)) {
      
      i = i + 1
      
      i.pu = pmgamma(u, mgshape, mgscale, mgweight)
      if (phiu) {
        i.phiu = 1 - i.pu
      } else {
        i.phiu = mean(x > u)
      }
      i.phib = (1 - i.phiu) / i.pu

      i.du = dmgamma(u, mgshape, mgscale, mgweight)
      sigmau = i.phiu/(i.phib * i.du)

      # Expectation step
      Mshape = matrix(mgshape, nrow = n, ncol = M, byrow = TRUE)
      Mscale = matrix(mgscale, nrow = n, ncol = M, byrow = TRUE)
      
      # unweighted density from each component
      Mdens = sapply(1:M, FUN = function(i) dgammagpd(x, mgshape[[i]], mgscale[[i]], u, sigmau, xi))
      
      # Reapply relative scaling phiu and phib
      Mdens[x > u,] = Mdens[x > u,] * i.phiu / (1 - i.pu)
      Mdens[x <= u,] = Mdens[x <= u,] * i.phib
      
      Mprob = t(mgweight * t(Mdens))
      
      # Estimated (posterior) probability of individual data being from each component,
      # conditional on each components' parameters
      tau = Mprob / rowSums(Mprob)

      # Updated mixture probabilities
      mgweight = colMeans(tau)
      
      # Maximisation step
      # MLE of each components' parameters,
      # conditional on probability of coming from each each component from E-step
      
      fit = try(optim(par = c(mgshape, mgscale, xi), fn = nluEMmgammagpdcon,
        u = u, tau = tau, mgweight = mgweight, x = x, M = M, phiu = phiu,
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
      xi = fit$par[2*M + 1]
      
      # store results
      nllh[i + 1] = nlmgammagpdcon(c(mgshape, mgscale, mgweight[-M], u, xi), x, M, phiu)    
      EMresults[i,] = c(nllh[i], mgshape, mgscale, mgweight, xi, fit$value, fit$convergence)    
    }
    
  } else { # complete (non-separable) likelihood
    
    # Initial values
    mgshape = pvector[1:M]
    mgscale = pvector[(M + 1):(2*M)]
    mgweight = pvector[(2*M + 1):(3*M - 1)]
    mgweight = c(mgweight, 1 - sum(mgweight)) # add constrained weight
    u = pvector[3*M]
    xi = pvector[3*M + 1]
    
    # check ordering of means, required for identifiability
    mui = mgshape * mgscale
    
    if (any(diff(mui) <= 0))
      stop("initial parameter vector does not satisfy constraint of strictly increasing component means")
    
    # Store iteration results
    nllh = rep(NA, maxit + 1)
    nllh[1] = 0
    
    EMresults = as.data.frame(matrix(NA, nrow = maxit, ncol = 3*M + 3 + 2))
    names(EMresults) = c("nllh", paste("mgshape", 1:M, sep = ""), paste("mgscale", 1:M, sep = ""),
      paste("mgweight", 1:M, sep = ""), "u", "xi", "nllh.cond.weights", "conv.cond.weights")
    
    # approximate negative log-likelihood at initial step
    i = 1
    nllh[i + 1] = nlmgammagpdcon(pvector, x, M, phiu)  

    if (is.infinite(nllh[i + 1])) stop("initial parameter values are invalid")
    
    EMresults[i,] = c(nllh[i + 1], mgshape, mgscale, mgweight, u, xi, NA, NA)

    while ((abs(nllh[i + 1] - nllh[i]) > abstol) & (i < maxit)) {

      i = i + 1
      
      i.pu = pmgamma(u, mgshape, mgscale, mgweight)
      if (phiu) {
        i.phiu = 1 - i.pu
      } else {
        i.phiu = mean(x > u)
      }
      i.phib = (1 - i.phiu) / i.pu

      i.du = dmgamma(u, mgshape, mgscale, mgweight)
      sigmau = i.phiu/(i.phib * i.du)

      # Expectation step
      Mshape = matrix(mgshape, nrow = n, ncol = M, byrow = TRUE)
      Mscale = matrix(mgscale, nrow = n, ncol = M, byrow = TRUE)
      
      # unweighted density from each component
      Mdens = sapply(1:M, FUN = function(i) dgammagpd(x, mgshape[[i]], mgscale[[i]], u, sigmau, xi))

      # Reapply relative scaling phiu and phib
      Mdens[x > u,] = Mdens[x > u,] * i.phiu / (1 - i.pu)
      Mdens[x <= u,] = Mdens[x <= u,] * i.phib
      
      Mprob = t(mgweight * t(Mdens))

      # Estimated (posterior) probability of individual data being from each component,
      # conditional on each components' parameters
      tau = Mprob / rowSums(Mprob)

      # Updated mixture probabilities
      mgweight = colMeans(tau)

      # Maximisation step
      # MLE of each components' parameters,
      # conditional on probability of coming from each each component from E-step
      
      fit = try(optim(par = c(mgshape, mgscale, u, xi), fn = nlEMmgammagpdcon,
        tau = tau, mgweight = mgweight, x = x, M = M, phiu = phiu,
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
      u = fit$par[2*M + 1]
      xi = fit$par[2*M + 2]

      # store results
      nllh[i + 1] = nlmgammagpdcon(c(mgshape, mgscale, mgweight[-M], u, xi), x, M, phiu)    
      EMresults[i,] = c(nllh[i], mgshape, mgscale, mgweight, u, xi, fit$value, fit$convergence)
    }
  }
    
  if (i == maxit) {
    warning("Maximum number of iteration reached, could be that M is too large or poor starting values")
  } else {
    print(paste(i, "iterations of EM algorithm"))
  }

  # Only output actual results, not blank rows
  EMresults = EMresults[1:i,]
  
  if ((fit$convergence != 0) | any(pvector[1:(3*M - 1)] == c(mgshape, mgscale, mgweight[-M])) |
      (pvector[length(pvector)] == xi) | (abs(fit$value) >= 1e6)) {
    conv = FALSE
    warning("check convergence")
  }

  if (!fixedu & (pvector[length(pvector) - 2] == u)) {
    conv = FALSE
    warning("check convergence")
  }

  pu = pmgamma(u, mgshape, mgscale, mgweight)
  if (phiu) {
    phiu = 1 - pu
    se.phiu = NA
  } else {
    phiu = mean(x > u, na.rm = TRUE)
    se.phiu = sqrt(phiu * (1 - phiu) / n)
  }
  phib = (1 - phiu) / pu

  du = dmgamma(u, mgshape, mgscale, mgweight)
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
  
  list(call = call, x = as.vector(x), init = as.vector(pvector), fixedu = fixedu, useq = useq, nllhuseq = nllhu,
    optim = fit, conv = conv, cov = invhess, mle = fit$par, se = se, rate = phiu,
    nllh = nllh[i + 1], n = n, M = M, mgshape = mgshape, mgscale = mgscale, mgweight = mgweight,
    u = u, sigmau = sigmau, xi = xi, phiu = phiu, se.phiu = se.phiu, EMresults = EMresults, posterior = tau)
}

#' @export
#' @aliases fmgammagpdcon lmgammagpdcon nlmgammagpdcon nlumgammagpdcon proflumgammagpdcon nlEMmgammagpdcon nluEMmgammagpdcon
#' @rdname  fmgammagpdcon

# log-likelihood function for mixture of gammas bulk with GPD for upper tail with continuity at threshold
lmgammagpdcon <- function(x, mgshape, mgscale, mgweight, u, xi, phiu = TRUE, log = TRUE) {
  
  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(u)
  check.param(xi)
  check.phiu(phiu, allowfalse = TRUE)
  check.logic(log)

  # user may try to input lists for mixture of gammas parameter vectors
  if (is.list(mgshape)) mgshape = unlist(mgshape)
  if (is.list(mgscale)) mgscale = unlist(mgscale)
  if (is.list(mgweight)) mgweight = unlist(mgweight)

  check.param(mgshape, allowvec = TRUE)
  check.param(mgscale, allowvec = TRUE)
  check.param(mgweight, allowvec = TRUE)

  # How many components in mixture
  M = check.inputn(c(length(mgshape), length(mgscale), length(mgweight)), allowscalar = TRUE)

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
  
  xu = x[which(x > u)]
  nu = length(xu)
  xb = x[which(x <= u)]
  nb = length(xb)
  if (n != nb + nu) {
    stop("total non-finite sample size is not equal to those above threshold and those below or equal to it")
  }

  # gamma components must be ordered to ensure identifiability (means strictly ascending order)
  mui = mgshape * mgscale

  if (any(mgscale <= 0) | any(mgshape <= 0) | any(diff(mui) <= 0) |
      any(mgweight <= 0) | (sum(mgweight) > 1) | (u <= 0) | (u <= min(x)) | (u >= max(x))) {
    l = -Inf
  } else {
    pu = pmgamma(u, mgshape, mgscale, mgweight)
    if (is.logical(phiu)) {
      if (phiu) {
        phiu = 1 - pu
      } else {
        phiu = nu/n
      }
    }
    phib = (1 - phiu)/pu

    du = dmgamma(u, mgshape, mgscale, mgweight)
    sigmau = phiu / (phib * du)
      
    syu = 1 + xi * (xu - u)/sigmau

    if ((min(syu) <= 0) | (sigmau <= 0) | (du < .Machine$double.eps) | (phiu <= 0) | (phiu >= 1) | (pu <= 0) | (pu >= 1)) {
      l = -Inf
    }
    else {
      l = lgpd(xu, u, sigmau, xi, phiu)
      l = l + sum(dmgamma(xb, mgshape, mgscale, mgweight, log = TRUE)) + nb * log(phib)
    }
  }
  
  if (!log) l = exp(l)
  
  l
}

#' @export
#' @aliases fmgammagpdcon lmgammagpdcon nlmgammagpdcon nlumgammagpdcon proflumgammagpdcon nlEMmgammagpdcon nluEMmgammagpdcon
#' @rdname  fmgammagpdcon

# negative log-likelihood function for mixture of gammas bulk with GPD for upper tail with continuity at threshold
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlmgammagpdcon <- function(pvector, x, M, phiu = TRUE, finitelik = FALSE) {

  check.n(M)
  if (M == 1) stop("use nlgammagpdcon instead")
  np = 3*M + 2 # maximum number of parameters

  # Check properties of inputs
  check.param(pvector, allowvec = TRUE)
  check.nparam(pvector, nparam = np - 1)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.phiu(phiu, allowfalse = TRUE)
  check.logic(finitelik)

  mgshape = pvector[1:M]
  mgscale = pvector[(M + 1):(2*M)]
  mgweight = pvector[(2*M + 1):(3*M - 1)]
  mgweight = c(mgweight, 1 - sum(mgweight)) # add constrained weight    
  u = pvector[3*M]
  xi = pvector[3*M + 1]
  
  nllh = -lmgammagpdcon(x, mgshape, mgscale, mgweight, u, xi, phiu) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}

#' @export
#' @aliases fmgammagpdcon lmgammagpdcon nlmgammagpdcon nlumgammagpdcon proflumgammagpdcon nlEMmgammagpdcon nluEMmgammagpdcon
#' @rdname  fmgammagpdcon

# negative log-likelihood function for mixture of gammas bulk with GPD for upper tail with continuity at threshold
# (wrapper for likelihood, designed for threshold to be fixed and other parameters optimised)
nlumgammagpdcon <- function(pvector, u, x, M, phiu = TRUE, finitelik = FALSE) {

  check.n(M)
  if (M == 1) stop("use nlugammagpdcon instead")
  np = 3*M + 2 # maximum number of parameters

  # Check properties of inputs
  check.param(pvector, allowvec = TRUE)
  check.nparam(pvector, nparam = np - 2)
  check.posparam(u)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.phiu(phiu, allowfalse = TRUE)
  check.logic(finitelik)

  mgshape = pvector[1:M]
  mgscale = pvector[(M + 1):(2*M)]
  mgweight = pvector[(2*M + 1):(3*M - 1)]
  mgweight = c(mgweight, 1 - sum(mgweight)) # add constrained weight    
  xi = pvector[3*M]
  
  nllh = -lmgammagpdcon(x, mgshape, mgscale, mgweight, u, xi, phiu) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}

#' @export
#' @aliases fmgammagpdcon lmgammagpdcon nlmgammagpdcon nlumgammagpdcon proflumgammagpdcon nlEMmgammagpdcon nluEMmgammagpdcon
#' @rdname  fmgammagpdcon

# negative log-likelihood function for mixture of gammas bulk with GPD for upper tail with continuity at threshold,
# with component probabilities separated for EM algorithm
nlEMmgammagpdcon <- function(pvector, tau, mgweight, x, M, phiu = TRUE, finitelik = FALSE) {

  check.n(M)
  if (M == 1) stop("use nlgammagpdcon instead")
  np.noweight = 2*M + 2 # maximum number of parameters

  # Check properties of inputs
  check.param(pvector, allowvec = TRUE)
  check.nparam(pvector, nparam = np.noweight)
  check.nparam(mgweight, nparam = M)
  check.prob(mgweight)
  
  if (abs(sum(mgweight) - 1) > 1e-6) stop("component weights must sum to one")

  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.phiu(phiu, allowfalse = TRUE)
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
  u = pvector[2*M + 1]
  xi = pvector[2*M + 2]
  
  # gamma components must be ordered to ensure identifiability (means strictly ascending order)
  mui = mgshape * mgscale

  if (any(mgscale <= 0) | any(mgshape <= 0) | any(diff(mui) <= 0) |
      (u <= 0) | (u <= min(x)) | (u >= max(x))) {
    nllh = Inf
  } else {
    pu = pmgamma(u, mgshape, mgscale, mgweight)
    if (is.logical(phiu)) {
      if (phiu) {
        phiu = 1 - pu
      } else {
        phiu = mean(x > u)
      }
    }
    phib = (1 - phiu) / pu
    
    du = dmgamma(u, mgshape, mgscale, mgweight)
    sigmau = phiu / (phib * du)
    
    Mshape = matrix(mgshape, nrow = n, ncol = M, byrow = TRUE)
    Mscale = matrix(mgscale, nrow = n, ncol = M, byrow = TRUE)

    # unweighted density from each component
    Mdens = sapply(1:M, FUN = function(i) dgammagpd(x, mgshape[[i]], mgscale[[i]], u, sigmau, xi))
    
    # Reapply relative scaling phiu and phib
    Mdens[x > u,] = Mdens[x > u,] * phiu / (1 - pu)
    Mdens[x <= u,] = Mdens[x <= u,] * phib
    
    nllh = -sum(tau * log(t(mgweight * t(Mdens))))
  }
    
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}

#' @export
#' @aliases fmgammagpdcon lmgammagpdcon nlmgammagpdcon nlumgammagpdcon proflumgammagpdcon nlEMmgammagpdcon nluEMmgammagpdcon
#' @rdname  fmgammagpdcon

# profile negative log-likelihood function for given threshold for
# mixture of gammas bulk with GPD for upper tail with continuity at threshold, using EM algorithm
# designed for sapply to loop over vector of thresholds (hence u is first input)
proflumgammagpdcon <- function(u, pvector, x, M, phiu = TRUE,
  method = "BFGS", control = list(maxit = 10000), finitelik = TRUE, ...) {

  check.n(M)
  if (M == 1) stop("use proflugammagpdcon instead")
  np = 3*M + 2 # maximum number of parameters

  # Check properties of inputs
  check.param(pvector, allowvec = TRUE)
  check.nparam(pvector, nparam = np - 2, allownull = TRUE)
  check.posparam(u)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.phiu(phiu, allowfalse = TRUE)
  check.logic(finitelik)
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

  if (is.unsorted(x)) {
    x = sort(x)
  } else {
    if (x[1] > x[length(x)])
      x = rev(x)
  }

  # approximate MLE for each gamma component
  approxmle = function(x) {
    s = log(mean(x)) - mean(log(x))
    k = (3 - s + sqrt((s - 3) ^ 2 + 24 * s))/12/s
    c(k, mean(x)/k)
  }

  if (!is.null(pvector)) {
    nllh = nlumgammagpdcon(pvector, u, x, M, phiu)
    if (is.infinite(nllh)) pvector = NULL
  }
  
  if (is.null(pvector)) {
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
    
    mgparam = simplify2array(tapply(x, INDEX = rep(1:M, times = nM), FUN = approxmle, simplify=TRUE))

    pvector = c(mgparam[1,], mgparam[2,], nM[-M]/n)
    
    # GPD defaults
    initfgpd = fgpd(x, u, std.err = FALSE)
    pvector[3*M] = initfgpd$xi
  }
  
  nllh = nlumgammagpdcon(pvector, u, x, M, phiu)

  # if still invalid then output cleanly
  if (is.infinite(nllh)) {
    warning(paste("initial parameter values for threshold u =", u, "are invalid"))
    fit = list(par = rep(NA, np), value = Inf, counts = 0, convergence = NA, 
      message = "initial values invalid", hessian = rep(NA, np))
  } else {

    # EM algorithm convergence conditions and constraints
    maxit = 1000
    abstol = 1e-8
    
    # Note each component is treated as gamma+gpd, where gpd is same in each but scaling
    # is phiu * mgweight and each gamma is scaled by phib * mgweight
    
    # Initial values
    mgshape = pvector[1:M]
    mgscale = pvector[(M + 1):(2*M)]
    mgweight = pvector[(2*M + 1):(3*M - 1)]
    mgweight = c(mgweight, 1 - sum(mgweight)) # add constrained weight
    xi = pvector[3*M]
    
    # check ordering of means, required for identifiability
    mui = mgshape * mgscale
    
    if (any(diff(mui) <= 0))
      stop("initial parameter vector does not satisfy constraint of strictly increasing component means")
    
    # Store iteration results
    nllh = rep(NA, maxit + 1)
    nllh[1] = 0
    
    # approximate negative log-likelihood at initial step
    i = 1
    nllh[i + 1] = nlumgammagpdcon(pvector, u, x, M, phiu)  
    
    if (is.infinite(nllh[i + 1])) stop("initial parameter values are invalid")
        
    while ((abs(nllh[i + 1] - nllh[i]) > abstol) & (i < maxit)) {
      
      i = i + 1
      
      i.pu = pmgamma(u, mgshape, mgscale, mgweight)
      if (is.logical(phiu)) {
        if (phiu) {
          i.phiu = 1 - i.pu
        } else {
          i.phiu = mean(x > u)
        }
      }
      i.phib = (1 - i.phiu) / i.pu

      i.du = dmgamma(u, mgshape, mgscale, mgweight)
      sigmau = i.phiu/(i.phib * i.du)

      # Expectation step
      Mshape = matrix(mgshape, nrow = n, ncol = M, byrow = TRUE)
      Mscale = matrix(mgscale, nrow = n, ncol = M, byrow = TRUE)
      
      # unweighted density from each component
      Mdens = sapply(1:M, FUN = function(i) dgammagpd(x, mgshape[[i]], mgscale[[i]], u, sigmau, xi))
      
      # Reapply relative scaling phiu and phib
      Mdens[x > u,] = Mdens[x > u,] * i.phiu / (1 - i.pu)
      Mdens[x <= u,] = Mdens[x <= u,] * i.phib
      
      Mprob = t(mgweight * t(Mdens))
      
      # Estimated (posterior) probability of individual data being from each component,
      # conditional on each components' parameters
      tau = Mprob / rowSums(Mprob)
      
      # Updated mixture probabilities
      mgweight = colMeans(tau)
      
      # Maximisation step
      # MLE of each components' parameters,
      # conditional on probability of coming from each each component from E-step
      
      fit = try(optim(par = c(mgshape, mgscale, xi), fn = nluEMmgammagpdcon,
        u = u, tau = tau, mgweight = mgweight, x = x, M = M, phiu = phiu,
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
      xi = fit$par[2*M + 1]
      
      # store results
      nllh[i + 1] = nlmgammagpdcon(c(mgshape, mgscale, mgweight[-M], u, xi), x, M, phiu)
    }
  }
    
  if (finitelik & is.infinite(nllh[i + 1])) {
    nllh[i + 1] = sign(nllh[i + 1]) * 1e6
  }

  nllh[i + 1]
}

#' @export
#' @aliases fmgammagpdcon lmgammagpdcon nlmgammagpdcon nlumgammagpdcon proflumgammagpdcon nlEMmgammagpdcon nluEMmgammagpdcon
#' @rdname  fmgammagpdcon

# negative log-likelihood function for mixture of gammas bulk with GPD for upper tail with continuity at threshold,
# with component probabilities separated for EM algorithm
# (wrapper for likelihood, designed for threshold to be fixed and other parameters optimised)
nluEMmgammagpdcon <- function(pvector, u, tau, mgweight, x, M, phiu = TRUE, finitelik = FALSE) {

  check.n(M)
  if (M == 1) stop("use nlugammagpdcon instead")
  np.noweight = 2*M + 2 # maximum number of parameters

  # Check properties of inputs
  check.param(pvector, allowvec = TRUE)
  check.nparam(pvector, nparam = np.noweight - 1)
  check.posparam(u)
  check.nparam(mgweight, nparam = M)
  check.prob(mgweight)
  
  if (abs(sum(mgweight) - 1) > 1e-6) stop("component weights must sum to one")
    
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.phiu(phiu, allowfalse = TRUE)
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
  xi = pvector[2*M + 1]
  
  # gamma components must be ordered to ensure identifiability (means strictly ascending order)
  mui = mgshape * mgscale

  if (any(mgscale <= 0) | any(mgshape <= 0) | any(diff(mui) <= 0) |
      (u <= 0) | (u <= min(x)) | (u >= max(x))) {
    nllh = Inf
  } else {
    pu = pmgamma(u, mgshape, mgscale, mgweight)
    if (is.logical(phiu)) {
      if (phiu) {
        phiu = 1 - pu
      } else {
        phiu = mean(x > u)
      }
    }
    phib = (1 - phiu)/pu

    du = dmgamma(u, mgshape, mgscale, mgweight)
    sigmau = phiu / (phib * du)
    
    Mshape = matrix(mgshape, nrow = n, ncol = M, byrow = TRUE)
    Mscale = matrix(mgscale, nrow = n, ncol = M, byrow = TRUE)

    # unweighted density from each component
    Mdens = sapply(1:M, FUN = function(i) dgammagpd(x, mgshape[[i]], mgscale[[i]], u, sigmau, xi))
    
    # Reapply relative scaling phiu and phib
    Mdens[x > u,] = Mdens[x > u,] * phiu / (1 - pu)
    Mdens[x <= u,] = Mdens[x <= u,] * phib
    
    nllh = -sum(tau * log(t(mgweight * t(Mdens))))
  }
    
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}
