#' @name kdengpd
#' 
#' @title Kernel Density Estimate and GPD Tail Extreme Value Mixture Model
#'
#' @description Density, cumulative distribution function, quantile function and
#'   random number generation for the extreme value mixture model with kernel density estimate for bulk
#'   distribution upto the threshold and conditional GPD above threshold. The parameters
#'   are the bandwidth \code{lambda}, threshold \code{u}
#'   GPD scale \code{sigmau} and shape \code{xi} and tail fraction \code{phiu}.
#'
#' @inheritParams normgpd
#' @inheritParams kden
#' @inheritParams gpd
#' 
#' @details Extreme value mixture model combining kernel density estimate (KDE) for the bulk
#' below the threshold and GPD for upper tail.
#' 
#' The user can pre-specify \code{phiu} 
#' permitting a parameterised value for the tail fraction \eqn{\phi_u}. Alternatively, when
#' \code{phiu=TRUE} the tail fraction is estimated as the tail fraction from the
#' KDE bulk model.
#' 
#' The alternate bandwidth definitions are discussed in the
#' \code{\link[evmix:kernels]{kernels}}, with the \code{lambda} as the default.
#' The \code{bw} specification is the same as used in the
#' \code{\link[stats:density]{density}} function.
#' 
#' The possible kernels are also defined in \code{\link[evmix:kernels]{kernels}}
#' with the \code{"gaussian"} as the default choice.
#' 
#' The cumulative distribution function with tail fraction \eqn{\phi_u} defined by the
#' upper tail fraction of the kernel density estimate (\code{phiu=TRUE}), upto the 
#' threshold \eqn{x \le u}, given by:
#' \deqn{F(x) = H(x)}
#' and above the threshold \eqn{x > u}:
#' \deqn{F(x) = H(u) + [1 - H(u)] G(x)}
#' where \eqn{H(x)} and \eqn{G(X)} are the KDE and conditional GPD
#' cumulative distribution functions respectively.
#' 
#' The cumulative distribution function for pre-specified \eqn{\phi_u}, upto the
#' threshold \eqn{x \le u}, is given by:
#' \deqn{F(x) = (1 - \phi_u) H(x)/H(u)}
#' and above the threshold \eqn{x > u}:
#' \deqn{F(x) = \phi_u + [1 - \phi_u] G(x)}
#' Notice that these definitions are equivalent when \eqn{\phi_u = 1 - H(u)}.
#' 
#' If no bandwidth is provided \code{lambda=NULL} and \code{bw=NULL} then the normal
#' reference rule is used, using the \code{\link[stats:bandwidth]{bw.nrd0}} function, which is
#' consistent with the \code{\link[stats:density]{density}} function. At least two kernel
#' centres must be provided as the variance needs to be estimated.
#' 
#' See \code{\link[evmix:gpd]{gpd}} for details of GPD upper tail component and 
#'\code{\link[evmix:kden]{dkden}} for details of KDE bulk component.
#' 
#' @return \code{\link[evmix:kdengpd]{dkdengpd}} gives the density, 
#' \code{\link[evmix:kdengpd]{pkdengpd}} gives the cumulative distribution function,
#' \code{\link[evmix:kdengpd]{qkdengpd}} gives the quantile function and 
#' \code{\link[evmix:kdengpd]{rkdengpd}} gives a random sample.
#' 
#' @note Unlike most of the other extreme value mixture model functions the 
#' \code{\link[evmix:kdengpd]{kdengpd}} functions have not been vectorised as
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
#' \code{\link[evmix:kdengpd]{rkdengpd}} is 1.
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
#' @aliases kdengpd dkdengpd pkdengpd qkdengpd rkdengpd
#' @family  kden kdengpd kdengpdcon bckden bckdengpd bckdengpdcon
#'          fkden fkdengpd fkdengpdcon fbckden fbckdengpd fbckdengpdcon
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' par(mfrow = c(2, 2))
#' 
#' kerncentres=rnorm(500, 0, 1)
#' xx = seq(-4, 4, 0.01)
#' hist(kerncentres, breaks = 100, freq = FALSE)
#' lines(xx, dkdengpd(xx, kerncentres, u = 1.2, sigmau = 0.56, xi = 0.1))
#' 
#' plot(xx, pkdengpd(xx, kerncentres), type = "l")
#' lines(xx, pkdengpd(xx, kerncentres, xi = 0.3), col = "red")
#' lines(xx, pkdengpd(xx, kerncentres, xi = -0.3), col = "blue")
#' legend("topleft", paste("xi =",c(0, 0.3, -0.3)),
#'       col=c("black", "red", "blue"), lty = 1, cex = 0.5)
#'
#' x = rkdengpd(1000, kerncentres, phiu = 0.1, u = 1.2, sigmau = 0.56, xi = 0.1)
#' xx = seq(-4, 6, 0.01)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-4, 6))
#' lines(xx, dkdengpd(xx, kerncentres, phiu = 0.1, u = 1.2, sigmau = 0.56, xi = 0.1))
#'
#' plot(xx, dkdengpd(xx, kerncentres, xi=0, phiu = 0.1), type = "l")
#' lines(xx, dkdengpd(xx, kerncentres, xi=0.2, phiu = 0.1), col = "red")
#' lines(xx, dkdengpd(xx, kerncentres, xi=-0.2, phiu = 0.1), col = "blue")
#' legend("topleft", c("xi = 0", "xi = 0.2", "xi = -0.2"),
#'       col=c("black", "red", "blue"), lty = 1)
#' }
#' 
NULL

#' @export
#' @aliases kdengpd dkdengpd pkdengpd qkdengpd rkdengpd
#' @rdname  kdengpd

# probability density function for kernel density estimate for the bulk
# distribution upto the threshold and conditional GPD above threshold.
dkdengpd <- function(x, kerncentres, lambda = NULL, u = as.vector(quantile(kerncentres, 0.9)), 
  sigmau = sqrt(6*var(kerncentres))/pi, xi = 0, phiu = TRUE,
  bw = NULL, kernel = "gaussian", log = FALSE) {
  
  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.quant(kerncentres, allowna = TRUE, allowinf = TRUE)
  check.kbw(lambda, bw, allownull = TRUE)
  check.kernel(kernel)
  check.param(u)
  check.posparam(sigmau)
  check.param(xi)
  check.phiu(phiu)
  check.logic(log)

  kernel = ifelse(kernel == "rectangular", "uniform", kernel)
  kernel = ifelse(kernel == "normal", "gaussian", kernel)

  if (any(is.infinite(x))) warning("infinite quantiles set to NA")

  x[is.infinite(x)] = NA # user will have to deal with infinite cases
    
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
  
  check.inputn(c(length(lambda), length(u), length(sigmau), length(xi), length(phiu)), allowscalar = TRUE) # scalar only
     
  pu = pkdenx(u, kerncentres, lambda, kernel)
  if (is.logical(phiu)) {
    phiu = 1 - pu
  } else {
    phiu = phiu
  }
  phib = (1 - phiu) / pu

  d = x # pass through NA/NaN as entered

  whichb = which(x <= u)
  nb = length(whichb)
  whichu = which(x > u)
  nu = length(whichu)

  if (nb > 0) d[whichb] = log(phib) + log(sapply(x[whichb], FUN = kdenx, kerncentres = kerncentres, lambda = lambda, kernel = kernel))
  if (nu > 0) d[whichu] = log(phiu) + dgpd(x[whichu], u, sigmau, xi, log = TRUE)
  
  if (!log) d = exp(d)
  
  d
}

#' @export
#' @aliases kdengpd dkdengpd pkdengpd qkdengpd rkdengpd
#' @rdname  kdengpd

# cumulative distribution function for kernel density estimate for the bulk
# distribution upto the threshold and conditional GPD above threshold.
pkdengpd <- function(q, kerncentres, lambda = NULL, u = as.vector(quantile(kerncentres, 0.9)), 
  sigmau = sqrt(6*var(kerncentres))/pi, xi = 0, phiu = TRUE,
  bw = NULL, kernel = "gaussian", lower.tail = TRUE) {
  
  # Check properties of inputs
  check.quant(q, allowna = TRUE, allowinf = TRUE)
  check.quant(kerncentres, allowna = TRUE, allowinf = TRUE)
  check.kbw(lambda, bw, allownull = TRUE)
  check.kernel(kernel)
  check.param(u)
  check.posparam(sigmau)
  check.param(xi)
  check.phiu(phiu)
  check.logic(lower.tail)

  kernel = ifelse(kernel == "rectangular", "uniform", kernel)
  kernel = ifelse(kernel == "normal", "gaussian", kernel)

  if (any(is.infinite(q))) warning("infinite quantiles set to NA")

  q[is.infinite(q)] = NA # user will have to deal with infinite cases
    
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
  
  check.inputn(c(length(lambda), length(u), length(sigmau), length(xi), length(phiu)), allowscalar = TRUE) # scalar only

  pu = pkdenx(u, kerncentres, lambda, kernel)
  if (is.logical(phiu)) {
    phiu = 1 - pu
  } else {
    phiu = phiu
  }
  phib = (1 - phiu) / pu
    
  p = q # pass through NA/NaN as entered
  
  whichb = which(q <= u)
  nb = length(whichb)
  whichu = which(q > u)
  nu = length(whichu)
  
  if (nb > 0) p[whichb] = phib*sapply(q[whichb], FUN = pkdenx, kerncentres = kerncentres, lambda = lambda, kernel = kernel)
  if (nu > 0) p[whichu] = (1 - phiu) + phiu*pgpd(q[whichu], u, sigmau, xi)
  
  if (!lower.tail) p = 1 - p
  
  p
}

#' @export
#' @aliases kdengpd dkdengpd pkdengpd qkdengpd rkdengpd
#' @rdname  kdengpd

# inverse cumulative distribution function for kernel density estimate for the bulk
# distribution upto the threshold and conditional GPD above threshold.
qkdengpd <- function(p, kerncentres, lambda = NULL, u = as.vector(quantile(kerncentres, 0.9)), 
  sigmau = sqrt(6*var(kerncentres))/pi, xi = 0, phiu = TRUE,
  bw = NULL, kernel = "gaussian", lower.tail = TRUE) {
  
  # Check properties of inputs
  check.prob(p, allowna = TRUE)
  check.quant(kerncentres, allowna = TRUE, allowinf = TRUE)
  check.kbw(lambda, bw, allownull = TRUE)
  check.kernel(kernel)
  check.param(u)
  check.posparam(sigmau)
  check.param(xi)
  check.phiu(phiu)
  check.logic(lower.tail)
    
  kernel = ifelse(kernel == "rectangular", "uniform", kernel)
  kernel = ifelse(kernel == "normal", "gaussian", kernel)

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
  
  check.inputn(c(length(lambda), length(u), length(sigmau), length(xi), length(phiu)), allowscalar = TRUE) # scalar only
    
  if (!lower.tail) p = 1 - p
  
  pu = pkdenx(u, kerncentres, lambda, kernel)
  if (is.logical(phiu)) {
    phiu = 1 - pu
  } else {
    phiu = phiu
  }
  phib = (1 - phiu) / pu
  
  q = p # pass through NA/NaN as entered

  whichb = which(p <= (1 - phiu))
  nb = length(whichb)
  whichu = which(p > (1 - phiu))
  nu = length(whichu)
  
  if (nb > 0) q[whichb] = qkden(p[whichb] / phib, kerncentres = kerncentres, lambda = lambda, kernel = kernel)
  if (nu > 0) q[whichu] = qgpd(p[whichu], u, sigmau, xi, phiu)
   
  q
}

#' @export
#' @aliases kdengpd dkdengpd pkdengpd qkdengpd rkdengpd
#' @rdname  kdengpd

# random number generation for kernel density estimate for the bulk
# distribution upto the threshold and conditional GPD above threshold.
rkdengpd <- function(n = 1, kerncentres, lambda = NULL, u = as.vector(quantile(kerncentres, 0.9)), 
  sigmau = sqrt(6*var(kerncentres))/pi, xi = 0, phiu = TRUE,
  bw = NULL, kernel = "gaussian") {
  
  # Check properties of inputs
  check.n(n)
  check.quant(kerncentres, allowna = TRUE, allowinf = TRUE)
  check.kbw(lambda, bw, allownull = TRUE)
  check.kernel(kernel)
  check.param(u)
  check.posparam(sigmau)
  check.param(xi)
  check.phiu(phiu)

  if (any(xi == 1)) stop("shape cannot be 1")
  
  qkdengpd(runif(n), kerncentres, lambda, u, sigmau, xi, phiu, bw, kernel)
}
