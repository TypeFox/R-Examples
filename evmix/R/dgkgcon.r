#' @name gkgcon
#' 
#' @title Kernel Density Estimate and GPD Both Upper and Lower Tails Extreme Value Mixture Model
#'  With Single Continuity Constraint at Both
#'
#' @description Density, cumulative distribution function, quantile function and
#'   random number generation for the extreme value mixture model with
#'   kernel density estimate for bulk distribution between thresholds and
#'   conditional GPD beyond thresholds and continuity at both of them. The parameters are the kernel bandwidth
#'  \code{lambda}, lower tail (threshold \code{ul}, 
#'   GPD shape \code{xil} and tail fraction \code{phiul})
#'   and upper tail (threshold \code{ur}, GPD shape 
#'   \code{xiR} and tail fraction \code{phiur}).
#'
#' @inheritParams gkg
#' 
#' @details Extreme value mixture model combining kernel density estimate (KDE) for the bulk
#' between thresholds and GPD beyond thresholds and continuity at both of them.
#' 
#' The user can pre-specify \code{phiul} and \code{phiur} 
#' permitting a parameterised value for the tail fractions \eqn{\phi_ul} and  \eqn{\phi_ur}.
#' Alternatively, when
#' \code{phiul=TRUE} and \code{phiur=TRUE} the tail fractions are estimated as the tail
#' fractions from the KDE bulk model.
#' 
#' The alternate bandwidth definitions are discussed in the
#' \code{\link[evmix:kernels]{kernels}}, with the \code{lambda} as the default.
#' The \code{bw} specification is the same as used in the
#' \code{\link[stats:density]{density}} function.
#' 
#' The possible kernels are also defined in \code{\link[evmix:kernels]{kernels}}
#' with the \code{"gaussian"} as the default choice.
#' 
#' Notice that the tail fraction cannot be 0 or 1, and the sum of upper and lower tail
#' fractions \code{phiul + phiur < 1}, so the lower threshold must be less than the upper, 
#' \code{ul < ur}.
#' 
#' The cumulative distribution function has three components. The lower tail with 
#' tail fraction \eqn{\phi_{ul}} defined by the KDE bulk model (\code{phiul=TRUE})
#' upto the lower threshold \eqn{x < u_l}:
#'   \deqn{F(x) = H(u_l) [1 - G_l(x)].}
#' where \eqn{H(x)} is the kernel density estimator cumulative distribution function (i.e. 
#' \code{mean(pnorm(x, kerncentres, bw))} and  
#' \eqn{G_l(X)} is the conditional GPD cumulative distribution function with negated
#' \eqn{x} value and threshold, i.e. \code{pgpd(-x, -ul, sigmaul, xil, phiul)}. The KDE
#' bulk model between the thresholds \eqn{u_l \le x \le u_r} given by:
#' \deqn{F(x) = H(x).}
#' Above the threshold \eqn{x > u_r} the usual conditional GPD:
#' \deqn{F(x) = H(u_r) + [1 - H(u_r)] G_r(x)}
#' where \eqn{G_r(X)} is the GPD cumulative distribution function, 
#' i.e. \code{pgpd(x, ur, sigmaur, xir, phiur)}.
#' 
#' The cumulative distribution function for the pre-specified tail fractions 
#' \eqn{\phi_{ul}} and \eqn{\phi_{ur}} is more complicated.  The unconditional GPD
#' is used for the lower tail \eqn{x < u_l}:
#'   \deqn{F(x) = \phi_{ul} [1 - G_l(x)].}
#' The KDE bulk model between the thresholds \eqn{u_l \le x \le u_r} given by:
#' \deqn{F(x) = \phi_{ul}+ (1-\phi_{ul}-\phi_{ur}) (H(x) - H(u_l)) / (H(u_r) - H(u_l)).}
#' Above the threshold \eqn{x > u_r} the usual conditional GPD:
#' \deqn{F(x) = (1-\phi_{ur}) + \phi_{ur} G(x)}
#' Notice that these definitions are equivalent when \eqn{\phi_{ul} = H(u_l)} and
#' \eqn{\phi_{ur} = 1 - H(u_r)}.
#' 
#' The continuity constraint at \code{ur} means that:
#'  \deqn{\phi_{ur} g_r(x) = (1-\phi_{ul}-\phi_{ur}) h(u_l)/ (H(u_r) - H(u_l)).}
#' By rearrangement, the GPD scale parameter \code{sigmaur} is then:
#' \deqn{\sigma_ur = \phi_{ur} (H(u_r) - H(u_l))/ h(u_l) (1-\phi_{ul}-\phi_{ur}).}
#' where \eqn{h(x)}, \eqn{g_l(x)} and \eqn{g_r(x)} are the KDE and conditional GPD
#' density functions for lower and upper tail respectively. 
#' In the special case of where the tail fraction is defined by the bulk model this reduces to
#' \deqn{\sigma_ur = [1-H(u_r)] / h(u_r)}.
#' 
#' The continuity constraint at \code{ul} means that:
#'  \deqn{\phi_{ul} g_l(x) = (1-\phi_{ul}-\phi_{ur}) h(u_l)/ (H(u_r) - H(u_l)).}
#' The GPD scale parameter \code{sigmaul} is replaced by:
#' \deqn{\sigma_ul = \phi_{ul} (H(u_r) - H(u_l))/ h(u_l) (1-\phi_{ul}-\phi_{ur}).}
#' In the special case of where the tail fraction is defined by the bulk model this reduces to
#' \deqn{\sigma_ul = H(u_l)/ h(u_l)}. 
#' 
#' If no bandwidth is provided \code{lambda=NULL} and \code{bw=NULL} then the normal
#' reference rule is used, using the \code{\link[stats:bandwidth]{bw.nrd0}} function, which is
#' consistent with the \code{\link[stats:density]{density}} function. At least two kernel
#' centres must be provided as the variance needs to be estimated.
#' 
#' See \code{\link[evmix:gpd]{gpd}} for details of GPD upper tail component and 
#'\code{\link[evmix:kden]{dkden}} for details of KDE bulk component.
#' 
#' @return \code{\link[evmix:gkgcon]{dgkgcon}} gives the density, 
#' \code{\link[evmix:gkgcon]{pgkgcon}} gives the cumulative distribution function,
#' \code{\link[evmix:gkgcon]{qgkgcon}} gives the quantile function and 
#' \code{\link[evmix:gkgcon]{rgkgcon}} gives a random sample.
#' 
#' @note Unlike most of the other extreme value mixture model functions the 
#' \code{\link[evmix:gkgcon]{gkgcon}} functions have not been vectorised as
#' this is not appropriate. The main inputs (\code{x}, \code{p} or \code{q})
#' must be either a scalar or a vector, which also define the output length.
#' The \code{kerncentres} can also be a scalar or vector.
#' 
#' The kernel centres \code{kerncentres} can either be a single datapoint or a vector
#' of data. The kernel centres (\code{kerncentres}) and locations to evaluate density (\code{x})
#' and cumulative distribution function (\code{q}) would usually be different.
#' 
#' Default values are provided for all inputs, except for the fundamentals 
#' \code{kerncentres}, \code{x}, \code{q} and \code{p}. The default sample size for 
#' \code{\link[evmix:gkgcon]{rgkgcon}} is 1.
#' 
#' Missing (\code{NA}) and Not-a-Number (\code{NaN}) values in \code{x},
#' \code{p} and \code{q} are passed through as is and infinite values are set to
#' \code{NA}. None of these are not permitted for the parameters or kernel centres.
#' 
#' Due to symmetry, the lower tail can be described by GPD by negating the quantiles. 
#' 
#' Error checking of the inputs (e.g. invalid probabilities) is carried out and
#' will either stop or give warning message as appropriate.
#' 
#' @references
#' \url{http://en.wikipedia.org/wiki/Kernel_density_estimation}
#' 
#' \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
#' 
#' Scarrott, C.J. and MacDonald, A. (2012). A review of extreme value
#' threshold estimation and uncertainty quantification. REVSTAT - Statistical
#' Journal 10(1), 33-59. Available from \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}
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
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}.
#'
#' @section Acknowledgments: Based on code
#' by Anna MacDonald produced for MATLAB.
#' 
#' @seealso \code{\link[evmix:kernels]{kernels}}, \code{\link[evmix:kfun]{kfun}},
#' \code{\link[stats:density]{density}}, \code{\link[stats:bandwidth]{bw.nrd0}}
#' and \code{\link[ks:kde.1d]{dkde}} in \code{\link[ks:kde.1d]{ks}} package.
#' 
#' @aliases gkgcon dgkgcon pgkgcon qgkgcon rgkgcon
#' @family  kden kdengpd kdengpdcon gkg gkgcon bckden bckdengpd bckdengpdcon
#'          fkden fkdengpd fkdengpdcon fgkg fgkgcon fbckden fbckdengpd fbckdengpdcon
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' par(mfrow = c(2, 2))
#' 
#' kerncentres=rnorm(1000,0,1)
#' x = rgkgcon(1000, kerncentres, phiul = 0.15, phiur = 0.15)
#' xx = seq(-6, 6, 0.01)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-6, 6))
#' lines(xx, dgkgcon(xx, kerncentres, phiul = 0.15, phiur = 0.15))
#' 
#' # three tail behaviours
#' plot(xx, pgkgcon(xx, kerncentres), type = "l")
#' lines(xx, pgkgcon(xx, kerncentres,xil = 0.3, xir = 0.3), col = "red")
#' lines(xx, pgkgcon(xx, kerncentres,xil = -0.3, xir = -0.3), col = "blue")
#' legend("topleft", paste("Symmetric xil=xir=",c(0, 0.3, -0.3)),
#'   col=c("black", "red", "blue"), lty = 1)
#' 
#' # asymmetric tail behaviours
#' x = rgkgcon(1000, kerncentres, xil = -0.3, phiul = 0.1, xir = 0.3, phiur = 0.1)
#' xx = seq(-6, 6, 0.01)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-6, 6))
#' lines(xx, dgkgcon(xx, kerncentres, xil = -0.3, phiul = 0.1, xir = 0.3, phiur = 0.1))
#' 
#' plot(xx, dgkgcon(xx, kerncentres, xil = -0.3, phiul = 0.2, xir = 0.3, phiur = 0.2),
#'   type = "l", ylim = c(0, 0.4))
#' lines(xx, dgkgcon(xx, kerncentres, xil = -0.3, phiul = 0.3, xir = 0.3, phiur = 0.3),
#'   col = "red")
#' lines(xx, dgkgcon(xx, kerncentres, xil = -0.3, phiul = TRUE, xir = 0.3, phiur = TRUE),
#'   col = "blue")
#' legend("topleft", c("phiul = phiur = 0.2", "phiul = phiur = 0.3", "Bulk Tail Fraction"),
#'   col=c("black", "red", "blue"), lty = 1)
#' }
#' 
NULL

#' @export
#' @aliases gkgcon dgkgcon pgkgcon qgkgcon rgkgcon
#' @rdname  gkgcon

# probability density function for kernel density estimate for the bulk
# distribution upto the thresholds and conditional GPD beyond thresholds
# with continuity at both
dgkgcon <- function(x, kerncentres, lambda = NULL,
  ul = as.vector(quantile(kerncentres, 0.1)), xil = 0, phiul = TRUE,
  ur = as.vector(quantile(kerncentres, 0.9)), xir = 0, phiur = TRUE,
  bw = NULL, kernel = "gaussian", log = FALSE) {
  
  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.quant(kerncentres, allowna = TRUE, allowinf = TRUE)
  check.kbw(lambda, bw, allownull = TRUE)
  check.kernel(kernel)
  check.param(ul)
  check.param(xil)
  check.phiu(phiul)
  check.param(ur)
  check.param(xir)
  check.phiu(phiur)
  check.logic(log)

  kernel = ifelse(kernel == "rectangular", "uniform", kernel)
  kernel = ifelse(kernel == "normal", "gaussian", kernel)

  if (any(is.infinite(x))) warning("infinite quantiles set to NA")

  x[is.infinite(x)] = NA # user will have to deal with infinite cases
    
  if (any(ul >= ur)) stop("lower threshold must be below upper threshold")
  
  if (!is.logical(phiul) & !is.logical(phiur)) {
    if ((phiul + phiur) > 1) stop("phiu + phiur must be less than 1")
  }

  if (any(!is.finite(kerncentres))) warning("non-finite kernel centres are dropped")

  kerncentres = kerncentres[is.finite(kerncentres)]
  check.quant(kerncentres)
  nk = length(kerncentres)

  if (is.null(lambda) & is.null(bw)) {
    if (nk == 1) {
      stop("Automated bandwidth estimation requires 2 or more kernel centres")
    } else if (nk < 10) {
      warning("Automated bandwidth estimation unreliable with less than 10 kernel centres")
    }
    bw = bw.nrd0(kerncentres)
  }
  lambda = klambda(bw, kernel, lambda)
  
  check.inputn(c(length(lambda), length(ul), length(xil), length(ur), length(xir),
    length(phiul), length(phiur)), allowscalar = TRUE) # scalar only
     
  if (is.logical(phiul)) {
    phiul = pkdenx(ul, kerncentres, lambda, kernel)
  } else {
    phiul = phiul
  }
  if (is.logical(phiur)) {
    phiur = 1 - pkdenx(ur, kerncentres, lambda, kernel)
  } else {
    phiur = phiur
  }
  phib = (1 - phiul - phiur) /
    (pkdenx(ur, kerncentres, lambda, kernel) - pkdenx(ul, kerncentres, lambda, kernel))

  sigmaul = phiul / (phib * kdenx(ul, kerncentres, lambda, kernel))
  sigmaur = phiur / (phib * kdenx(ur, kerncentres, lambda, kernel)) 
  
  check.posparam(sigmaul)
  check.posparam(sigmaur)
    
  dgkg(x, kerncentres, lambda, ul, sigmaul, xil, phiul, ur, sigmaur, xir, phiur,
    kernel = kernel, log = log)
}

#' @export
#' @aliases gkgcon dgkgcon pgkgcon qgkgcon rgkgcon
#' @rdname  gkgcon

# cumulative distribution function for kernel density estimate for the bulk
# distribution upto the thresholds and conditional GPD beyond thresholds
# with continuity at both
pgkgcon <- function(q, kerncentres, lambda = NULL,
  ul = as.vector(quantile(kerncentres, 0.1)), xil = 0, phiul = TRUE,
  ur = as.vector(quantile(kerncentres, 0.9)), xir = 0, phiur = TRUE,
  bw = NULL, kernel = "gaussian", lower.tail = TRUE) {
  
  # Check properties of inputs
  check.quant(q, allowna = TRUE, allowinf = TRUE)
  check.quant(kerncentres, allowna = TRUE, allowinf = TRUE)
  check.kbw(lambda, bw, allownull = TRUE)
  check.kernel(kernel)
  check.param(ul)
  check.param(xil)
  check.phiu(phiul)
  check.param(ur)
  check.param(xir)
  check.phiu(phiur)
  check.logic(lower.tail)

  kernel = ifelse(kernel == "rectangular", "uniform", kernel)
  kernel = ifelse(kernel == "normal", "gaussian", kernel)

  if (any(is.infinite(q))) warning("infinite quantiles set to NA")

  q[is.infinite(q)] = NA # user will have to deal with infinite cases
    
  if (any(ul >= ur)) stop("lower threshold must be below upper threshold")
  
  if (!is.logical(phiul) & !is.logical(phiur)) {
    if ((phiul + phiur) > 1) stop("phiu + phiur must be less than 1")
  }

  if (any(!is.finite(kerncentres))) warning("non-finite kernel centres are dropped")

  kerncentres = kerncentres[is.finite(kerncentres)]
  check.quant(kerncentres)
  nk = length(kerncentres)

  if (is.null(lambda) & is.null(bw)) {
    if (nk == 1) {
      stop("Automated bandwidth estimation requires 2 or more kernel centres")
    } else if (nk < 10) {
      warning("Automated bandwidth estimation unreliable with less than 10 kernel centres")
    }
    bw = bw.nrd0(kerncentres)
  }
  lambda = klambda(bw, kernel, lambda)
  
  check.inputn(c(length(lambda), length(ul), length(xil), length(ur), length(xir),
    length(phiul), length(phiur)), allowscalar = TRUE) # scalar only
     
  if (is.logical(phiul)) {
    phiul = pkdenx(ul, kerncentres, lambda, kernel)
  } else {
    phiul = phiul
  }
  if (is.logical(phiur)) {
    phiur = 1 - pkdenx(ur, kerncentres, lambda, kernel)
  } else {
    phiur = phiur
  }
  phib = (1 - phiul - phiur) /
    (pkdenx(ur, kerncentres, lambda, kernel) - pkdenx(ul, kerncentres, lambda, kernel))
    
  sigmaul = phiul / (phib * kdenx(ul, kerncentres, lambda, kernel))
  sigmaur = phiur / (phib * kdenx(ur, kerncentres, lambda, kernel)) 
  
  check.posparam(sigmaul)
  check.posparam(sigmaur)
    
  pgkg(q, kerncentres, lambda, ul, sigmaul, xil, phiul, ur, sigmaur, xir, phiur,
    kernel = kernel, lower.tail = lower.tail)
}

#' @export
#' @aliases gkgcon dgkgcon pgkgcon qgkgcon rgkgcon
#' @rdname  gkgcon

# inverse cumulative distribution function for kernel density estimate for the bulk
# distribution upto the thresholds and conditional GPD beyond thresholds
# with continuity at both
qgkgcon <- function(p, kerncentres, lambda = NULL,
  ul = as.vector(quantile(kerncentres, 0.1)), xil = 0, phiul = TRUE,
  ur = as.vector(quantile(kerncentres, 0.9)), xir = 0, phiur = TRUE,
  bw = NULL, kernel = "gaussian", lower.tail = TRUE) {
  
  # Check properties of inputs
  check.prob(p, allowna = TRUE)
  check.quant(kerncentres, allowna = TRUE, allowinf = TRUE)
  check.kbw(lambda, bw, allownull = TRUE)
  check.kernel(kernel)
  check.param(ul)
  check.param(xil)
  check.phiu(phiul)
  check.param(ur)
  check.param(xir)
  check.phiu(phiur)
  check.logic(lower.tail)
    
  kernel = ifelse(kernel == "rectangular", "uniform", kernel)
  kernel = ifelse(kernel == "normal", "gaussian", kernel)

  if (any(ul >= ur)) stop("lower threshold must be below upper threshold")
  
  if (!is.logical(phiul) & !is.logical(phiur)) {
    if ((phiul + phiur) > 1) stop("phiu + phiur must be less than 1")
  }

  if (any(!is.finite(kerncentres))) warning("non-finite kernel centres are dropped")

  kerncentres = kerncentres[is.finite(kerncentres)]
  check.quant(kerncentres)
  nk = length(kerncentres)

  if (is.null(lambda) & is.null(bw)) {
    if (nk == 1) {
      stop("Automated bandwidth estimation requires 2 or more kernel centres")
    } else if (nk < 10) {
      warning("Automated bandwidth estimation unreliable with less than 10 kernel centres")
    }
    bw = bw.nrd0(kerncentres)
  }
  lambda = klambda(bw, kernel, lambda)
  
  check.inputn(c(length(lambda), length(ul), length(xil), length(ur), length(xir),
    length(phiul), length(phiur)), allowscalar = TRUE) # scalar only
     
  if (is.logical(phiul)) {
    phiul = pkdenx(ul, kerncentres, lambda, kernel)
  } else {
    phiul = phiul
  }
  if (is.logical(phiur)) {
    phiur = 1 - pkdenx(ur, kerncentres, lambda, kernel)
  } else {
    phiur = phiur
  }
  phib = (1 - phiul - phiur) /
    (pkdenx(ur, kerncentres, lambda, kernel) - pkdenx(ul, kerncentres, lambda, kernel))
  
  sigmaul = phiul / (phib * kdenx(ul, kerncentres, lambda, kernel))
  sigmaur = phiur / (phib * kdenx(ur, kerncentres, lambda, kernel)) 
  
  check.posparam(sigmaul)
  check.posparam(sigmaur)
    
  qgkg(p, kerncentres, lambda, ul, sigmaul, xil, phiul, ur, sigmaur, xir, phiur,
    kernel = kernel, lower.tail = lower.tail)
}

#' @export
#' @aliases gkgcon dgkgcon pgkgcon qgkgcon rgkgcon
#' @rdname  gkgcon

# random number generation for kernel density estimate for the bulk
# distribution upto the thresholds and conditional GPD beyond thresholds
# with continuity at both
rgkgcon <- function(n = 1, kerncentres, lambda = NULL,
  ul = as.vector(quantile(kerncentres, 0.1)), xil = 0, phiul = TRUE,
  ur = as.vector(quantile(kerncentres, 0.9)), xir = 0, phiur = TRUE,
  bw = NULL, kernel = "gaussian") {
  
  # Check properties of inputs
  check.n(n)
  check.quant(kerncentres, allowna = TRUE, allowinf = TRUE)
  check.kbw(lambda, bw, allownull = TRUE)
  check.kernel(kernel)
  check.param(ul)
  check.param(xil)
  check.phiu(phiul)
  check.param(ur)
  check.param(xir)
  check.phiu(phiur)
     
  if (any(xil == 1) | any(xir == 1)) stop("shape cannot be 1")

  qgkgcon(runif(n), kerncentres, lambda, ul, xil, phiul, ur, xir, phiur, bw, kernel)
}
