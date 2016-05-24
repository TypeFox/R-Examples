#' @name gkg
#' 
#' @title Kernel Density Estimate and GPD Both Upper and Lower Tails Extreme Value Mixture Model
#'
#' @description Density, cumulative distribution function, quantile function and
#'   random number generation for the extreme value mixture model with kernel density estimate for bulk
#'   distribution between thresholds and conditional GPD beyond thresholds. The parameters are the kernel bandwidth
#'  \code{lambda}, lower tail (threshold \code{ul}, 
#'   GPD scale \code{sigmaul} and shape \code{xil} and tail fraction \code{phiul})
#'   and upper tail (threshold \code{ur}, GPD scale \code{sigmaur} and shape 
#'   \code{xiR} and tail fraction \code{phiur}).
#'
#' @inheritParams gng
#' @inheritParams kdengpd
#' @inheritParams kden
#' @inheritParams kernels
#' @inheritParams gpd
#' 
#' @details Extreme value mixture model combining kernel density estimate (KDE) for the bulk
#' between thresholds and GPD beyond thresholds.
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
#' If no bandwidth is provided \code{lambda=NULL} and \code{bw=NULL} then the normal
#' reference rule is used, using the \code{\link[stats:bandwidth]{bw.nrd0}} function, which is
#' consistent with the \code{\link[stats:density]{density}} function. At least two kernel
#' centres must be provided as the variance needs to be estimated.
#' 
#' See \code{\link[evmix:gpd]{gpd}} for details of GPD upper tail component and 
#'\code{\link[evmix:kden]{dkden}} for details of KDE bulk component.
#' 
#' @return \code{\link[evmix:gkg]{dgkg}} gives the density, 
#' \code{\link[evmix:gkg]{pgkg}} gives the cumulative distribution function,
#' \code{\link[evmix:gkg]{qgkg}} gives the quantile function and 
#' \code{\link[evmix:gkg]{rgkg}} gives a random sample.
#' 
#' @note Unlike most of the other extreme value mixture model functions the 
#' \code{\link[evmix:gkg]{gkg}} functions have not been vectorised as
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
#' \code{\link[evmix:gkg]{rgkg}} is 1.
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
#' @aliases gkg dgkg pgkg qgkg rgkg
#' @family  kden kdengpd kdengpdcon gkg gkgcon bckden bckdengpd bckdengpdcon
#'          fkden fkdengpd fkdengpdcon fgkg fgkgcon fbckden fbckdengpd fbckdengpdcon
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' par(mfrow = c(2, 2))
#' 
#' kerncentres=rnorm(1000,0,1)
#' x = rgkg(1000, kerncentres, phiul = 0.15, phiur = 0.15)
#' xx = seq(-6, 6, 0.01)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-6, 6))
#' lines(xx, dgkg(xx, kerncentres, phiul = 0.15, phiur = 0.15))
#' 
#' # three tail behaviours
#' plot(xx, pgkg(xx, kerncentres), type = "l")
#' lines(xx, pgkg(xx, kerncentres,xil = 0.3, xir = 0.3), col = "red")
#' lines(xx, pgkg(xx, kerncentres,xil = -0.3, xir = -0.3), col = "blue")
#' legend("topleft", paste("Symmetric xil=xir=",c(0, 0.3, -0.3)),
#'   col=c("black", "red", "blue"), lty = 1)
#' 
#' # asymmetric tail behaviours
#' x = rgkg(1000, kerncentres, xil = -0.3, phiul = 0.1, xir = 0.3, phiur = 0.1)
#' xx = seq(-6, 6, 0.01)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-6, 6))
#' lines(xx, dgkg(xx, kerncentres, xil = -0.3, phiul = 0.1, xir = 0.3, phiur = 0.1))
#' 
#' plot(xx, dgkg(xx, kerncentres, xil = -0.3, phiul = 0.2, xir = 0.3, phiur = 0.2),
#'   type = "l", ylim = c(0, 0.4))
#' lines(xx, dgkg(xx, kerncentres, xil = -0.3, phiul = 0.3, xir = 0.3, phiur = 0.3),
#'   col = "red")
#' lines(xx, dgkg(xx, kerncentres, xil = -0.3, phiul = TRUE, xir = 0.3, phiur = TRUE),
#'   col = "blue")
#' legend("topleft", c("phiul = phiur = 0.2", "phiul = phiur = 0.3", "Bulk Tail Fraction"),
#'   col=c("black", "red", "blue"), lty = 1)
#' }
#' 
NULL

#' @export
#' @aliases gkg dgkg pgkg qgkg rgkg
#' @rdname  gkg

# probability density function for kernel density estimate for the bulk
# distribution upto the thresholds and conditional GPD beyond thresholds
dgkg <- function(x, kerncentres, lambda = NULL,
  ul = as.vector(quantile(kerncentres, 0.1)), sigmaul = sqrt(6*var(kerncentres))/pi, xil = 0, phiul = TRUE,
  ur = as.vector(quantile(kerncentres, 0.9)), sigmaur = sqrt(6*var(kerncentres))/pi, xir = 0, phiur = TRUE,
  bw = NULL, kernel = "gaussian", log = FALSE) {
  
  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.quant(kerncentres, allowna = TRUE, allowinf = TRUE)
  check.kbw(lambda, bw, allownull = TRUE)
  check.kernel(kernel)
  check.param(ul)
  check.posparam(sigmaul)
  check.param(xil)
  check.phiu(phiul)
  check.param(ur)
  check.posparam(sigmaur)
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
  
  check.inputn(c(length(lambda),
    length(ul), length(sigmaul), length(xil),
    length(ur), length(sigmaur), length(xir), length(phiul), length(phiur)), allowscalar = TRUE) # scalar only
     
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

  d = x # pass through NA/NaN as entered

  whichul = which(x < ul)
  nul = length(whichul)
  whichb = which((x <= ur) & (x >= ul)) 
  nb = length(whichb)
  whichur = which(x > ur)
  nur = length(whichur)
    
  if (nul > 0) d[whichul] = log(phiul) + dgpd(-x[whichul], -ul, sigmaul, xil, log = TRUE)
  if (nb > 0) d[whichb] = log(phib) + log(sapply(x[whichb], FUN = kdenx, kerncentres = kerncentres, lambda = lambda, kernel = kernel))
  if (nur > 0) d[whichur] = log(phiur) + dgpd(x[whichur], ur, sigmaur, xir, log = TRUE)
  
  if (!log) d = exp(d)
  
  d
}

#' @export
#' @aliases gkg dgkg pgkg qgkg rgkg
#' @rdname  gkg

# cumulative distribution function for kernel density estimate for the bulk
# distribution upto the thresholds and conditional GPD beyond thresholds
pgkg <- function(q, kerncentres, lambda = NULL,
  ul = as.vector(quantile(kerncentres, 0.1)), sigmaul = sqrt(6*var(kerncentres))/pi, xil = 0, phiul = TRUE,
  ur = as.vector(quantile(kerncentres, 0.9)), sigmaur = sqrt(6*var(kerncentres))/pi, xir = 0, phiur = TRUE,
  bw = NULL, kernel = "gaussian", lower.tail = TRUE) {
  
  # Check properties of inputs
  check.quant(q, allowna = TRUE, allowinf = TRUE)
  check.quant(kerncentres, allowna = TRUE, allowinf = TRUE)
  check.kbw(lambda, bw, allownull = TRUE)
  check.kernel(kernel)
  check.param(ul)
  check.posparam(sigmaul)
  check.param(xil)
  check.phiu(phiul)
  check.param(ur)
  check.posparam(sigmaur)
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
  
  check.inputn(c(length(lambda),
    length(ul), length(sigmaul), length(xil),
    length(ur), length(sigmaur), length(xir), length(phiul), length(phiur)), allowscalar = TRUE) # scalar only
     
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
    
  p = q # pass through NA/NaN as entered
  
  whichul = which(q < ul)
  nul = length(whichul)
  whichb = which((q <= ur) & (q >= ul)) 
  nb = length(whichb)
  whichur = which(q > ur)
  nur = length(whichur)
  
  if (nul > 0) p[whichul] = 1 - pgpd(-q[whichul], -ul, sigmaul, xil, phiul)
  if (nb > 0) p[whichb] = phiul + phib*(sapply(q[whichb], FUN = pkdenx, kerncentres = kerncentres, lambda = lambda, kernel = kernel) - pkdenx(ul, kerncentres, lambda, kernel))
  if (nur > 0) p[whichur] = pgpd(q[whichur], ur, sigmaur, xir, phiur)
  
  if (!lower.tail) p = 1 - p
  
  p
}

#' @export
#' @aliases gkg dgkg pgkg qgkg rgkg
#' @rdname  gkg

# inverse cumulative distribution function for kernel density estimate for the bulk
# distribution upto the thresholds and conditional GPD beyond thresholds
qgkg <- function(p, kerncentres, lambda = NULL,
  ul = as.vector(quantile(kerncentres, 0.1)), sigmaul = sqrt(6*var(kerncentres))/pi, xil = 0, phiul = TRUE,
  ur = as.vector(quantile(kerncentres, 0.9)), sigmaur = sqrt(6*var(kerncentres))/pi, xir = 0, phiur = TRUE,
  bw = NULL, kernel = "gaussian", lower.tail = TRUE) {
  
  # Check properties of inputs
  check.prob(p, allowna = TRUE)
  check.quant(kerncentres, allowna = TRUE, allowinf = TRUE)
  check.kbw(lambda, bw, allownull = TRUE)
  check.kernel(kernel)
  check.param(ul)
  check.posparam(sigmaul)
  check.param(xil)
  check.phiu(phiul)
  check.param(ur)
  check.posparam(sigmaur)
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
  
  check.inputn(c(length(lambda),
    length(ul), length(sigmaul), length(xil),
    length(ur), length(sigmaur), length(xir), length(phiul), length(phiur)), allowscalar = TRUE) # scalar only
     
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
  
  q = p # pass through NA/NaN as entered

  whichul = which(p < phiul)
  nul = length(whichul)
  whichb = which((p <= (1 - phiur)) & (p >= phiul)) 
  nb = length(whichb)
  whichur = which(p > (1 - phiur))
  nur = length(whichur)
  
  if (nul > 0) q[whichul] = -qgpd(1 - p[whichul], -ul, sigmaul, xil, phiul)
  if (nb > 0) q[whichb] = qkden((p[whichb] - phiul) / phib + pkdenx(ul, kerncentres, lambda, kernel), kerncentres, lambda, kernel = kernel)
  if (nur > 0) q[whichur] = qgpd(p[whichur], ur, sigmaur, xir, phiur)
   
  q
}

#' @export
#' @aliases gkg dgkg pgkg qgkg rgkg
#' @rdname  gkg

# random number generation for kernel density estimate for the bulk
# distribution upto the thresholds and conditional GPD beyond thresholds
rgkg <- function(n = 1, kerncentres, lambda = NULL,
  ul = as.vector(quantile(kerncentres, 0.1)), sigmaul = sqrt(6*var(kerncentres))/pi, xil = 0, phiul = TRUE,
  ur = as.vector(quantile(kerncentres, 0.9)), sigmaur = sqrt(6*var(kerncentres))/pi, xir = 0, phiur = TRUE,
  bw = NULL, kernel = "gaussian") {
  
  # Check properties of inputs
  check.n(n)
  check.quant(kerncentres, allowna = TRUE, allowinf = TRUE)
  check.kbw(lambda, bw, allownull = TRUE)
  check.kernel(kernel)
  check.param(ul)
  check.posparam(sigmaul)
  check.param(xil)
  check.phiu(phiul)
  check.param(ur)
  check.posparam(sigmaur)
  check.param(xir)
  check.phiu(phiur)
     
  if (any(xil == 1) | any(xir == 1)) stop("shape cannot be 1")

  qgkg(runif(n), kerncentres, lambda, ul, sigmaul, xil, phiul, ur, sigmaur, xir, phiur, bw, kernel)
}
