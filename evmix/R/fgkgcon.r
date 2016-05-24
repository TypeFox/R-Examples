#' @export
#' 
#' @title MLE Fitting of Kernel Density Estimate for Bulk and GPD for Both Tails with 
#'  Single Continuity Constraint at Both Thresholds Extreme Value Mixture Model
#'
#' @description Maximum likelihood estimation for fitting the extreme value 
#' mixture model with kernel density estimate for bulk distribution between thresholds and conditional
#' GPDs for both tails with continuity at thresholds. With options for profile likelihood estimation for both thresholds and
#' fixed threshold approach.
#'
#' @inheritParams fgkg
#' 
#' @details The extreme value mixture model with kernel density estimate for bulk and
#' GPD for both tails with continuity at thresholds is 
#' fitted to the entire dataset using maximum likelihood estimation. The estimated
#' parameters, variance-covariance matrix and their standard errors are automatically
#' output.
#' 
#' See help for \code{\link[evmix:fnormgpd]{fnormgpd}} and \code{\link[evmix:fgkg]{fgng}} 
#' for details, type \code{help fnormgpd} and \code{help fgng}. 
#' Only the different features are outlined below for brevity.
#' 
#' The GPD \code{sigmaul} and \code{sigmaur} parameters are now specified as function of
#' other parameters, see 
#' help for \code{\link[evmix:gkgcon]{dgkgcon}} for details, type \code{help gkgcon}.
#' Therefore, \code{sigmaul} and \code{sigmaur} should not be included in the parameter
#' vector if initial values are provided, making the full parameter vector 
#' The full parameter vector is
#' (\code{lambda}, \code{ul}, \code{xil}, \code{ur}, \code{xir})
#' if thresholds are also estimated and
#' (\code{lambda}, \code{xil}, \code{xir})
#' for profile likelihood or fixed threshold approach.
#' 
#' Cross-validation likelihood is used for KDE, but standard likelihood is used
#' for GPD components. See help for \code{\link[evmix:fkden]{fkden}} for details,
#' type \code{help fkden}.
#' 
#' The alternate bandwidth definitions are discussed in the 
#' \code{\link[evmix:kernels]{kernels}}, with the \code{lambda} as the default
#' used in the likelihood fitting. The \code{bw} specification is the same as
#' used in the \code{\link[stats:density]{density}} function.
#' 
#' The possible kernels are also defined in \code{\link[evmix:kernels]{kernels}}
#' with the \code{"gaussian"} as the default choice.
#' 
#' The tail fractions \code{phiul} and \code{phiur} are treated separately to the other parameters, 
#' to allow for all their representations. In the fitting functions 
#' \code{\link[evmix:fgkgcon]{fgkgcon}} and
#' \code{\link[evmix:fgkgcon]{proflugkgcon}} they are logical:
#' \itemize{
#'  \item default values \code{phiul=TRUE} and \code{phiur=TRUE} - tail fractions specified by 
#'    KDE distribution and survivior functions respectively and
#'    standard error is output as \code{NA}.
#'  \item \code{phiul=FALSE} and \code{phiur=FALSE} - treated as extra parameters estimated using
#'    the MLE which is the sample proportion beyond the thresholds and 
#'    standard error is output.
#' }
#' In the likelihood functions \code{\link[evmix:fgkgcon]{lgkgcon}},
#' \code{\link[evmix:fgkgcon]{nlgkgcon}} and \code{\link[evmix:fgkgcon]{nlugkgcon}} 
#' it can be logical or numeric:
#' \itemize{
#'  \item logical - same as for fitting functions with default values \code{phiul=TRUE} and \code{phiur=TRUE}.
#'  \item numeric - any value over range \eqn{(0, 1)}. Notice that the tail
#'    fraction probability cannot be 0 or 1 otherwise there would be no
#'    contribution from either tail or bulk components respectively. Also,
#'    \code{phiul+phiur<1} as bulk must contribute.
#' }
#' 
#' If the profile likelihood approach is used, then a grid search over all combinations of both thresholds
#' is carried out. The combinations which lead to less than 5 in any datapoints beyond the thresholds are not considered.
#' 
#' @section Warning:
#' See important warnings about cross-validation likelihood estimation in 
#' \code{\link[evmix:fkden]{fkden}}, type \code{help fkden}.
#' 
#' @return Log-likelihood is given by \code{\link[evmix:fgkgcon]{lgkgcon}} and it's
#'   wrappers for negative log-likelihood from \code{\link[evmix:fgkgcon]{nlgkgcon}}
#'   and \code{\link[evmix:fgkgcon]{nlugkgcon}}. Profile likelihood for both
#'   thresholds given by \code{\link[evmix:fgkgcon]{proflugkgcon}}. Fitting function
#'   \code{\link[evmix:fgkgcon]{fgkgcon}} returns a simple list with the
#'   following elements
#'
#' \tabular{ll}{
#'  \code{call}:      \tab \code{optim} call\cr
#'  \code{x}:         \tab data vector \code{x}\cr
#'  \code{init}:      \tab \code{pvector}\cr
#'  \code{fixedu}:    \tab fixed thresholds, logical\cr
#'  \code{ulseq}:     \tab lower threshold vector for profile likelihood or scalar for fixed threshold\cr
#'  \code{urseq}:     \tab upper threshold vector for profile likelihood or scalar for fixed threshold\cr
#'  \code{nllhuseq}:  \tab profile negative log-likelihood at each threshold pair in (ulseq, urseq)\cr
#'  \code{optim}:     \tab complete \code{optim} output\cr
#'  \code{mle}:       \tab vector of MLE of parameters\cr
#'  \code{cov}:       \tab variance-covariance matrix of MLE of parameters\cr
#'  \code{se}:        \tab vector of standard errors of MLE of parameters\cr
#'  \code{rate}:      \tab \code{phiu} to be consistent with \code{\link[evd:fpot]{evd}}\cr
#'  \code{nllh}:      \tab minimum negative log-likelihood\cr
#'  \code{n}:         \tab total sample size\cr
#'  \code{lambda}:    \tab MLE of lambda (kernel half-width)\cr
#'  \code{ul}:        \tab lower threshold (fixed or MLE)\cr
#'  \code{sigmaul}:   \tab MLE of lower tail GPD scale (estimated from other parameters)\cr
#'  \code{xil}:       \tab MLE of lower tail GPD shape\cr
#'  \code{phiul}:     \tab MLE of lower tail fraction (bulk model or parameterised approach)\cr
#'  \code{se.phiul}:  \tab standard error of MLE of lower tail fraction\cr
#'  \code{ur}:        \tab upper threshold (fixed or MLE)\cr
#'  \code{sigmaur}:   \tab MLE of upper tail GPD scale (estimated from other parameters)\cr
#'  \code{xir}:       \tab MLE of upper tail GPD shape\cr
#'  \code{phiur}:     \tab MLE of upper tail fraction (bulk model or parameterised approach)\cr
#'  \code{se.phiur}:  \tab standard error of MLE of lower tail fraction\cr
#'  \code{bw}:        \tab MLE of bw (kernel standard deviations)\cr
#'  \code{kernel}:    \tab kernel name\cr
#' }
#' 
#' @note The data and kernel centres are both vectors. Infinite and missing sample values
#' (and kernel centres) are dropped.
#' 
#' When \code{pvector=NULL} then the initial values are:
#' \itemize{
#'  \item normal reference rule for bandwidth, using the \code{\link[stats:bandwidth]{bw.nrd0}} function, which is
#'    consistent with the \code{\link[stats:density]{density}} function. At least two kernel
#'    centres must be provided as the variance needs to be estimated.
#'  \item lower threshold 10\% quantile (not relevant for profile likelihood for threshold or fixed threshold approaches);
#'  \item upper threshold 90\% quantile (not relevant for profile likelihood for threshold or fixed threshold approaches);
#'  \item MLE of GPD shape parameters beyond thresholds. 
#' }
#' 
#' @references
#' \url{http://www.math.canterbury.ac.nz/~c.scarrott/evmix}
#' 
#' \url{http://en.wikipedia.org/wiki/Kernel_density_estimation}
#' 
#' \url{http://en.wikipedia.org/wiki/Cross-validation_(statistics)}
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
#' Bowman, A.W. (1984). An alternative method of cross-validation for the smoothing of
#' density estimates. Biometrika 71(2), 353-360.
#' 
#' Duin, R.P.W. (1976). On the choice of smoothing parameters for Parzen estimators of
#' probability density functions. IEEE Transactions on Computers C25(11), 1175-1179.
#' 
#' MacDonald, A., Scarrott, C.J., Lee, D., Darlow, B., Reale, M. and Russell, G. (2011).
#' A flexible extreme value mixture model. Computational Statistics and Data Analysis
#' 55(6), 2137-2157.
#' 
#' Wand, M. and Jones, M.C. (1995). Kernel Smoothing. Chapman && Hall.
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @section Acknowledgments: See Acknowledgments in
#'   \code{\link[evmix:fnormgpd]{fnormgpd}}, type \code{help fnormgpd}. Based on code
#' by Anna MacDonald produced for MATLAB.
#' 
#' @seealso \code{\link[evmix:kernels]{kernels}}, \code{\link[evmix:kfun]{kfun}},
#'  \code{\link[stats:density]{density}}, \code{\link[stats:bandwidth]{bw.nrd0}}
#' and \code{\link[ks:kde.1d]{dkde}} in \code{\link[ks:kde.1d]{ks}} package.
#'  \code{\link[evmix:fgpd]{fgpd}} and \code{\link[evmix:gpd]{gpd}}
#'  
#' @aliases fgkgcon lgkgcon nlgkgcon proflugkgcon nlugkgcon
#' @family  kdengpd kdengpdcon fkdengpd fkdengpdcon normgpd fnormgpd gkg gkgcon fgkg fgkgcon
#'          kden bckden bckdengpd bckdengpdcon fkden fbckden fbckdengpd fbckdengpdcon
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
#' fit = fgkgcon(x)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-4, 4))
#' lines(xx, y)
#' with(fit, lines(xx, dgkgcon(xx, x, lambda, ul, xil, phiul,
#'    ur, xir, phiur), col="red"))
#' abline(v = c(fit$ul, fit$ur), col = "red")
#'   
#' # No continuity constraint
#' fit2 = fgkg(x)
#' with(fit2, lines(xx, dgkg(xx, x, lambda, ul, sigmaul, xil, phiul,
#'    ur, sigmaur, xir, phiur), col="blue"))
#' abline(v = c(fit2$ul, fit2$ur), col = "blue")
#' legend("topleft", c("True Density","No continuity constraint","With continuty constraint"),
#'   col=c("black", "blue", "red"), lty = 1)
#'   
#' # Profile likelihood for initial value of threshold and fixed threshold approach
#' fitu = fgkgcon(x, ulseq = seq(-2, -0.2, length = 10), 
#'  urseq = seq(0.2, 2, length = 10))
#' fitfix = fgkgcon(x, ulseq = seq(-2, -0.2, length = 10), 
#'  urseq = seq(0.2, 2, length = 10), fixedu = TRUE)
#' 
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-4, 4))
#' lines(xx, y)
#' with(fit, lines(xx, dgkgcon(xx, x, lambda, ul, xil, phiul,
#'    ur, xir, phiur), col="red"))
#' abline(v = c(fit$ul, fit$ur), col = "red")
#' with(fitu, lines(xx, dgkgcon(xx, x, lambda, ul, xil, phiul,
#'    ur, xir, phiur), col="purple"))
#' abline(v = c(fitu$ul, fitu$ur), col = "purple")
#' with(fitfix, lines(xx, dgkgcon(xx, x, lambda, ul, xil, phiul,
#'    ur, xir, phiur), col="darkgreen"))
#' abline(v = c(fitfix$ul, fitfix$ur), col = "darkgreen")
#' legend("topright", c("True Density","Default initial value (90% quantile)",
#'  "Prof. lik. for initial value", "Prof. lik. for fixed threshold"),
#'  col=c("black", "red", "purple", "darkgreen"), lty = 1)
#' }
#'   

# maximum likelihood fitting for kernel density estimate for bulk with GPD for
# both tails with continuity at thresholds
fgkgcon <- function(x, phiul = TRUE, phiur = TRUE, ulseq = NULL, urseq = NULL, fixedu = FALSE,
  pvector = NULL, kernel = "gaussian", add.jitter = FALSE, factor = 0.1, amount = NULL,
  std.err = TRUE, method = "BFGS", control = list(maxit = 10000), finitelik = TRUE, ...) {

  call <- match.call()
    
  np = 5 # maximum number of parameters

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.logic(phiul)
  check.logic(phiur)
  check.param(ulseq, allowvec = TRUE, allownull = TRUE)
  check.param(urseq, allowvec = TRUE, allownull = TRUE)
  check.logic(fixedu)
  check.logic(std.err)
  check.optim(method)
  check.control(control)
  check.logic(finitelik)

  check.kernel(kernel)
  check.posparam(factor)
  check.posparam(amount, allownull = TRUE)
  check.logic(add.jitter)
  
  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    x = x[is.finite(x)] # ignore missing and infinite cases
  }

  check.quant(x)
  n = length(x)

  if (add.jitter) x = jitter(x, factor, amount)

  xuniq = unique(x)
  if (length(xuniq) < (0.95*n))
    warning("data may be rounded, as more than 5% are ties, so bandwidth could be biased to zero")

  if ((method == "L-BFGS-B") | (method == "BFGS")) finitelik = TRUE
  
  # useq must be specified if threshold is fixed
  if (fixedu & (is.null(ulseq) | is.null(urseq)))
    stop("for fixed threshold approach, ulseq and urseq must be specified (as scalar or vector)")
  
  # Check if profile likelihood or fixed threshold is being used
  # and determine initial values for parameters in each case
  if (is.null(ulseq) | is.null(ulseq)) { # not profile or fixed

    check.nparam(pvector, nparam = np, allownull = TRUE)
    
    if (is.null(pvector)) {
      if (n == 1) {
        stop("Automated bandwidth estimation requires 2 or more kernel centres")
      } else if (n < 10) {
        stop("Automated bandwidth estimation unreliable with less than 10 kernel centres")
      }
      pvector[1] = klambda(bw.nrd0(x), kernel)
      pvector[2] = as.vector(quantile(x, 0.1))
      initfgpd = fgpd(-x, -pvector[2], std.err = FALSE)
      pvector[3] = initfgpd$xi
      pvector[4] = as.vector(quantile(x, 0.9))
      initfgpd = fgpd(x, pvector[4], std.err = FALSE)
      pvector[5] = initfgpd$xi
    }
    
  } else { # profile or fixed
    
    check.nparam(pvector, nparam = np - 2, allownull = TRUE)

    # profile likelihood for threshold or scalar given
    if ((length(ulseq) != 1) | (length(urseq) != 1)) {
      
      # remove thresholds with less than 5 excesses
      ulseq = ulseq[sapply(ulseq, FUN = function(u, x) sum(x < u) > 5, x = x)]
      check.param(ulseq, allowvec = TRUE)
      urseq = urseq[sapply(urseq, FUN = function(u, x) sum(x > u) > 5, x = x)]
      check.param(urseq, allowvec = TRUE)

      ulrseq = expand.grid(ulseq, urseq)
      
      # remove those where ulseq >= urseq
      if (any(ulrseq[1] >= ulrseq[2])) {
        warning("lower thresholds above or equal to upper threshold are ignored")
        ulrseq = ulrseq[ulrseq[1] < ulrseq[2],]
      }

      nllhu = apply(ulrseq, 1, proflugkgcon, pvector = pvector, x = x,
        phiul = phiul, phiur = phiur, kernel = kernel,
        method = method, control = control, finitelik = finitelik, ...)
      
      if (all(!is.finite(nllhu))) stop("thresholds are all invalid")
      ul = ulrseq[which.min(nllhu), 1]
      ur = ulrseq[which.min(nllhu), 2]

    } else {
      if (ulseq >= urseq) stop("lower threshold cannot be above or equal to upper threshold")
      ul = ulseq
      ur = urseq
    }

    if (fixedu) { # threshold fixed
      if (is.null(pvector)) {
        if (n == 1) {
          stop("Automated bandwidth estimation requires 2 or more kernel centres")
        } else if (n < 10) {
          stop("Automated bandwidth estimation unreliable with less than 10 kernel centres")
        }
        pvector[1] = klambda(bw.nrd0(x), kernel)
        initfgpd = fgpd(-x, -ul, std.err = FALSE)
        pvector[2] = initfgpd$xi
        initfgpd = fgpd(x, ur, std.err = FALSE)
        pvector[3] = initfgpd$xi
      }
    } else { # threshold as initial value in usual MLE
      if (is.null(pvector)) {
        if (n == 1) {
          stop("Automated bandwidth estimation requires 2 or more kernel centres")
        } else if (n < 10) {
          stop("Automated bandwidth estimation unreliable with less than 10 kernel centres")
        }
        pvector[1] = klambda(bw.nrd0(x), kernel)
        pvector[2] = ul
        initfgpd = fgpd(-x, -pvector[2], std.err = FALSE)
        pvector[3] = initfgpd$xi
        pvector[4] = ur
        initfgpd = fgpd(x, pvector[4], std.err = FALSE)
        pvector[5] = initfgpd$xi
      } else {
        pvector[5] = pvector[3] # shift upper tail GPD shape to add in ur
        pvector[4] = ur
        pvector[3] = pvector[2] # shift lower tail GPD shape to add in ul
        pvector[2] = ul
      }
    }
  }

  if (fixedu) { # fixed threshold (separable) likelihood
    nllh = nlugkgcon(pvector, ul, ur, x, phiul, phiur, kernel = kernel)
    if (is.infinite(nllh)) {
      pvector[c(2, 3)] = 0.1
      nllh = nlugkgcon(pvector, ul, ur, x, phiul, phiur, kernel = kernel)    
    }
    if (is.infinite(nllh)) stop("initial parameter values are invalid")
  
    fit = optim(par = as.vector(pvector), fn = nlugkgcon, ul = ul, ur = ur, x = x,
      phiul = phiul, phiur = phiur, kernel = kernel,
      finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)    
    
    lambda = fit$par[1]
    xil = fit$par[2]
    xir = fit$par[3]
    
  } else { # complete (non-separable) likelihood
    
    nllh = nlgkgcon(pvector, x, phiul, phiur, kernel = kernel)
    if (is.infinite(nllh)) {
      pvector[c(3, 5)] = 0.1
      nllh = nlgkgcon(pvector, x, phiul, phiur, kernel = kernel)    
    }
    if (is.infinite(nllh)) stop("initial parameter values are invalid")
  
    fit = optim(par = as.vector(pvector), fn = nlgkgcon, x = x, phiul = phiul, phiur = phiur, kernel = kernel,
      finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)    
    
    lambda = fit$par[1]
    ul = fit$par[2]
    xil = fit$par[3]
    ur = fit$par[4]
    xir = fit$par[5]
  }

  bw = kbw(fit$par[1], kernel)

  conv = TRUE
  if ((fit$convergence != 0) | any(fit$par == pvector) | (abs(fit$value) >= 1e6)) {
    conv = FALSE
    warning("check convergence")
  }

  pul = pkdenx(ul, x, lambda, kernel)
  if (phiul) {
    phiul = pul
    se.phiul = NA
  } else {
    phiul = mean(x < ul, na.rm = TRUE)
    se.phiul = sqrt(phiul * (1 - phiul) / n)
  }
  pur = pkdenx(ur, x, lambda, kernel)
  if (phiur) {
    phiur = 1 - pur
    se.phiur = NA
  } else {
    phiur = mean(x > ur, na.rm = TRUE)
    se.phiur = sqrt(phiur * (1 - phiur) / n)
  }
  phib = (1 - phiul - phiur) / (pur - pul)
    
  dul = kdenx(ul, x, lambda, kernel)
  dur = kdenx(ur, x, lambda, kernel)
    
  sigmaul = phiul / (phib * dul)
  sigmaur = phiur / (phib * dur)
  
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
    init = as.vector(pvector), fixedu = fixedu, ulseq = ulseq, urseq = urseq, nllhuseq = nllhu,
    optim = fit, conv = conv, cov = invhess, mle = fit$par, se = se, ratel = phiul, rater = phiur,
    nllh = fit$value, n = n, lambda = lambda, 
    ul = ul, sigmaul = sigmaul, xil = xil, phiul = phiul, se.phiul = se.phiul, 
    ur = ur, sigmaur = sigmaur, xir = xir, phiur = phiur, se.phiur = se.phiur, bw = bw, kernel = kernel)
}

#' @export
#' @aliases fgkgcon lgkgcon nlgkgcon proflugkgcon nlugkgcon
#' @rdname  fgkgcon

# log-likelihood function for kernel density estimate for bulk with GPD for
# both tails with continuity at thresholds
# cross-validation for KDE component
lgkgcon <- function(x, lambda = NULL, ul = 0, xil = 0, phiul = TRUE,
  ur = 0, xir = 0, phiur = TRUE, bw = NULL, kernel = "gaussian", log = TRUE) {

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(lambda, allownull = TRUE)
  check.param(bw, allownull = TRUE)
  check.param(ul)                       
  check.param(xil)
  check.param(ur)                       
  check.param(xir)
  check.phiu(phiul, allowfalse = TRUE)
  check.phiu(phiur, allowfalse = TRUE)
  check.logic(log)

  check.kernel(kernel)

  kernel = ifelse(kernel == "rectangular", "uniform", kernel)
  kernel = ifelse(kernel == "normal", "gaussian", kernel)
  
  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    x = x[is.finite(x)] # ignore missing and infinite cases
  }

  check.quant(x)
  n = length(x)

  if (is.null(lambda) & is.null(bw)) {
    if (n == 1) {
      stop("Automated bandwidth estimation requires 2 or more kernel centres")
    } else if (n < 10) {
      warning("Automated bandwidth estimation unreliable with less than 10 kernel centres")
    }
    bw = bw.nrd0(x)
  }

  if (is.null(lambda)) lambda = klambda(bw, kernel, lambda)
  
  check.inputn(c(length(lambda), length(ul), length(xil), length(phiul), 
    length(ur), length(xir), length(phiur)), allowscalar = TRUE)

  # assume NA or NaN are irrelevant as entire lower tail is now modelled
  # inconsistent with evd library definition
  # hence use which() to ignore these

  xur = x[which(x > ur)]
  nur = length(xur)
  xul = x[which(x < ul)]
  nul = length(xul)
  xb = x[which((x >= ul) & (x <= ur))]
  nb = length(xb)

  if ((lambda <= 0) | (ul <= min(x)) | (ul >= max(x)) | (ur <= min(x)) | (ur >= max(x))| (ur <= ul)) {
    l = -Inf
  } else {
    if (is.logical(phiul)) {
      pul = pkdenx(ul, x, lambda, kernel)
      if (phiul) {
        phiul = pul
      } else {
        phiul = nul / n
      }
    }
    if (is.logical(phiur)) {
      pur = pkdenx(ur, x, lambda, kernel)
      if (phiur) {
        phiur = 1 - pur
      } else {
        phiur = nur / n
      }
    }
    phib = (1 - phiul - phiur) / (pur - pul)
  
    dul = kdenx(ul, x, lambda, kernel)
    dur = kdenx(ur, x, lambda, kernel)
    
    sigmaul = phiul / (phib * dul)
    sigmaur = phiur / (phib * dur) 

    syul = 1 + xil * (ul - xul) / sigmaul  
    syur = 1 + xir * (xur - ur) / sigmaur  
  
    if ((min(syul) <= 0) | (phiul <= 0) | (phiul >= 1) | 
        (min(syur) <= 0) | (phiur <= 0) | (phiur >= 1) | ((phiul + phiur) > 1) |
        (pul <= 0) | (pul >= 1) | (pur <= 0) | (pur >= 1) |
        (phib < .Machine$double.eps) |
        (dul < .Machine$double.eps) | (dur < .Machine$double.eps) | 
        (sigmaul <= 0) | (sigmaur <= 0)) {
      l = -Inf
    } else { 
      l = lgpd(-xul, -ul, sigmaul, xil, phiul)
      l = l + lgpd(xur, ur, sigmaur, xir, phiur)
      l = l + lkden(xb, lambda, kernel = kernel, extracentres = c(xul, xur)) + nb*log(phib)
    }
  }
  
  if (!log) l = exp(l)
  
  l
}

#' @export
#' @aliases fgkgcon lgkgcon nlgkgcon proflugkgcon nlugkgcon
#' @rdname  fgkgcon

# negative log-likelihood function for kernel density estimate for bulk with GPD for
# both tails with continuity at thresholds
# cross-validation for KDE component
# (wrapper for likelihood, inputs and checks designed for optimisation)
nlgkgcon <- function(pvector, x, phiul = TRUE, phiur = TRUE, kernel = "gaussian", finitelik = FALSE) {

  np = 5 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.phiu(phiul, allowfalse = TRUE)
  check.phiu(phiur, allowfalse = TRUE)
  check.logic(finitelik)

  check.kernel(kernel)
  kernel = ifelse(kernel == "rectangular", "uniform", kernel)
  kernel = ifelse(kernel == "normal", "gaussian", kernel)

  lambda = pvector[1]
  ul = pvector[2]
  xil = pvector[3]
  ur = pvector[4]
  xir = pvector[5]

  nllh = -lgkgcon(x, lambda, ul, xil, phiul, ur, xir, phiur, kernel = kernel) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}

#' @export
#' @aliases fgkgcon lgkgcon nlgkgcon proflugkgcon nlugkgcon
#' @rdname  fgkgcon

# profile negative log-likelihood function for given threshold for
# kernel density estimate for bulk with GPD for both tails with continuity at thresholds
# designed for apply to loop over vector of thresholds (hence c(ul, ur) vector is first input)
# cross-validation for KDE component
proflugkgcon <- function(ulr, pvector, x, phiul = TRUE, phiur = TRUE, kernel = "gaussian",
  method = "BFGS", control = list(maxit = 10000), finitelik = TRUE, ...) {

  np = 5 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np - 2, allownull = TRUE)
  check.param(ulr, allowvec = TRUE)
  check.nparam(ulr, nparam = 2)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.phiu(phiul, allowfalse = TRUE)
  check.phiu(phiur, allowfalse = TRUE)
  check.optim(method)
  check.control(control)
  check.logic(finitelik)

  check.kernel(kernel)

  kernel = ifelse(kernel == "rectangular", "uniform", kernel)
  kernel = ifelse(kernel == "normal", "gaussian", kernel)

  if (any(!is.finite(x))) {
    warning("non-finite cases have been removed")
    x = x[is.finite(x)] # ignore missing and infinite cases
  }

  check.quant(x)
  n = length(x)
  
  ul = ulr[1]
  ur = ulr[2]

  if (ul >= ur) stop("lower threshold cannot be above or equal to upper threshold")

  # check initial values for other parameters, try usual alternative
  if (!is.null(pvector)) {
    nllh = nlugkgcon(pvector, ul, ur, x, phiul, phiur, kernel = kernel)
    
    if (is.infinite(nllh)) pvector = NULL
  }

  if (is.null(pvector)) {
    if (n == 1) {
      stop("Automated bandwidth estimation requires 2 or more kernel centres")
    } else if (n < 10) {
        stop("Automated bandwidth estimation unreliable with less than 10 kernel centres")
    }
    pvector[1] = klambda(bw.nrd0(x), kernel)
    initfgpd = fgpd(-x, -ul, std.err = FALSE)
    pvector[2] = initfgpd$xi
    initfgpd = fgpd(x, ur, std.err = FALSE)
    pvector[3] = initfgpd$xi
    nllh = nlugkgcon(pvector, ul, ur, x, phiul, phiur, kernel = kernel)
  }

  if (is.infinite(nllh)) {
    pvector[c(2, 3)] = 0.1
    nllh = nlugkgcon(pvector, ul, ur, x, phiul, phiur, kernel = kernel)    
  }

  # if still invalid then output cleanly
  if (is.infinite(nllh)) {
    warning(paste("initial parameter values for thresholds ul =", ul, "and ur =", ur,"are invalid"))
    fit = list(par = rep(NA, np), value = Inf, counts = 0, convergence = NA, 
      message = "initial values invalid", hessian = rep(NA, np))
  } else {

    fit = optim(par = as.vector(pvector), fn = nlugkgcon, ul = ul, ur = ur, x = x, 
      phiul = phiul, phiur = phiur, kernel = kernel,
      finitelik = finitelik, method = method, control = control, hessian = TRUE, ...)
  }
    
  if (finitelik & is.infinite(fit$value)) {
    fit$value = sign(fit$value) * 1e6
  }

  fit$value
}

#' @export
#' @aliases fgkgcon lgkgcon nlgkgcon proflugkgcon nlugkgcon
#' @rdname  fgkgcon

# negative log-likelihood function for kernel density estimate for bulk with GPD for
# both tails with continuity at thresholds
# (wrapper for likelihood, designed for threshold to be fixed and other parameters optimised)
# cross-validation for KDE component
nlugkgcon <- function(pvector, ul, ur, x, phiul = TRUE, phiur = TRUE, kernel = "gaussian", finitelik = FALSE) {

  np = 5 # maximum number of parameters

  # Check properties of inputs
  check.nparam(pvector, nparam = np - 2)
  check.param(ul)
  check.param(ur)
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.phiu(phiul, allowfalse = TRUE)
  check.phiu(phiur, allowfalse = TRUE)
  check.logic(finitelik)

  check.kernel(kernel)
  kernel = ifelse(kernel == "rectangular", "uniform", kernel)
  kernel = ifelse(kernel == "normal", "gaussian", kernel)
    
  if (ul >= ur) stop("lower threshold cannot be above or equal to upper threshold")

  lambda = pvector[1]
  xil = pvector[2]
  xir = pvector[3]

  nllh = -lgkgcon(x, lambda, ul, xil, phiul, ur, xir, phiur, kernel = kernel) 
  
  if (finitelik & is.infinite(nllh)) {
    nllh = sign(nllh) * 1e6
  }

  nllh
}
