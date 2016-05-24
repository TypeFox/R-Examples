#' @name kden
#' 
#' @title Kernel Density Estimation, With Variety of Kernels
#'
#' @description Density, cumulative distribution function, quantile function and
#'   random number generation for the kernel density estimation using the kernel
#'   specified by \code{kernel}, with a constant bandwidth specified by either
#'   \code{lambda} or \code{bw}.
#'
#' @inheritParams gpd
#' @inheritParams kernels
#' 
#' @details Kernel density estimation using one of many possible kernels with a
#' constant bandwidth. 
#' 
#' The alternate bandwidth definitions are discussed in the
#' \code{\link[evmix:kernels]{kernels}}, with the \code{lambda} as the default.
#' The \code{bw} specification is the same as used in the
#' \code{\link[stats:density]{density}} function.
#' 
#' The possible kernels are also defined in \code{\link[evmix:kernels]{kernels}} help
#' documentation with the \code{"gaussian"} as the default choice.
#'
#' The density function \code{\link[evmix:kden]{dkden}} produces exactly the
#' same density estimate as \code{\link[stats:density]{density}} when a sequence
#' of \code{x} values are provided, see examples. The latter function is far
#' more efficient in this situation as it takes advantage of the computational
#' savings from doing the kernel smoothing in the spectral domain (using the FFT),
#' where the convolution becomes a multiplication. So even after accounting for applying
#' the (Fast) Fourier Transform (FFT) and its inverse it is much more efficient
#' especially for a large sample size or large number of evaluation points.
#' 
#' However, this KDE function applies the less efficient convolution using the
#' standard definition:
#' \deqn{\hat{f}_(x) = \frac{1}{n} \sum_{j=1}^{n} K(\frac{x - x_j}{\lambda})}
#' where \eqn{K(.)} is the density function for the standard
#' kernel. Thus are no restriction on the values \code{x} can take. For example, in the 
#' \code{"gaussian"} kernel case for a particular \code{x} the density is evaluated as
#' \code{mean(dnorm(x, kerncentres, lambda))} for the density and
#' \code{mean(pnorm(x, kerncentres, lambda))} for cumulative distribution
#' function which is slower than the FFT but is more adaptable.
#' 
#' An inversion sampler is used for random number generation which also rather
#' inefficient, as it can be carried out more efficiently using a mixture representation.
#' 
#' The quantile function is rather complicated as there is no closed form solution,
#' so is obtained by numerical approximation of the inverse cumulative distribution function
#' \eqn{P(X \le q) = p} to find \eqn{q}. The quantile function 
#' \code{\link[evmix:kden]{qkden}} evaluates the KDE cumulative distribution
#' function over the range from \code{c(max(kerncentre) - lambda, max(kerncentre) + lambda)},
#' or \code{c(max(kerncentre) - 5*lambda, max(kerncentre) + 5*lambda)} for normal kernel.
#' Outside of this range the quantiles are set to \code{-Inf} for lower tail and \code{Inf}
#' for upper tail. A sequence of values
#' of length fifty times the number of kernels (with minimum of 1000) is first
#' calculated. Spline based interpolation using \code{\link[stats:splinefun]{splinefun}},
#' with default \code{monoH.FC} method, is then used to approximate the quantile
#' function. This is a similar approach to that taken
#' by Matt Wand in the \code{\link[ks:kde.1d]{qkde}} in the \code{\link[ks:kde.1d]{ks}} package.
#' 
#' If no bandwidth is provided \code{lambda=NULL} and \code{bw=NULL} then the normal
#' reference rule is used, using the \code{\link[stats:bandwidth]{bw.nrd0}} function, which is
#' consistent with the \code{\link[stats:density]{density}} function. At least two kernel
#' centres must be provided as the variance needs to be estimated.
#' 
#' @return \code{\link[evmix:kden]{dkden}} gives the density, 
#' \code{\link[evmix:kden]{pkden}} gives the cumulative distribution function,
#' \code{\link[evmix:kden]{qkden}} gives the quantile function and 
#' \code{\link[evmix:kden]{rkden}} gives a random sample.
#' 
#' @note Unlike most of the other extreme value mixture model functions the 
#'   \code{\link[evmix:kden]{kden}} functions have not been vectorised as
#'   this is not appropriate. The main inputs (\code{x}, \code{p} or \code{q})
#'   must be either a scalar or a vector, which also define the output length.
#' 
#' The kernel centres \code{kerncentres} can either be a single datapoint or a vector
#' of data. The kernel centres (\code{kerncentres}) and locations to evaluate density (\code{x})
#' and cumulative distribution function (\code{q}) would usually be different.
#' 
#' Default values are provided for all inputs, except for the fundamentals 
#' \code{kerncentres}, \code{x}, \code{q} and \code{p}. The default sample size for 
#' \code{\link[evmix:kden]{rkden}} is 1.
#' 
#' Missing (\code{NA}) and Not-a-Number (\code{NaN}) values in \code{x},
#' \code{p} and \code{q} are passed through as is and infinite values are set to
#' \code{NA}. None of these are not permitted for the parameters.
#' 
#' Error checking of the inputs (e.g. invalid probabilities) is carried out and
#' will either stop or give warning message as appropriate.
#' 
#' @references
#' \url{http://en.wikipedia.org/wiki/Kernel_density_estimation}
#' 
#' \url{http://en.wikipedia.org/wiki/Cross-validation_(statistics)}
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
#' @aliases kden dkden pkden qkden rkden
#' @family  kden kdengpd kdengpdcon bckden bckdengpd bckdengpdcon
#'          fkden fkdengpd fkdengpdcon fbckden fbckdengpd fbckdengpdcon
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' par(mfrow = c(2, 2))
#' 
#' nk=50
#' x = rnorm(nk)
#' xx = seq(-5, 5, 0.01)
#' plot(xx, dnorm(xx))
#' rug(x)
#' for (i in 1:nk) lines(xx, dnorm(xx, x[i], sd = bw.nrd0(x))*0.05)
#' lines(xx, dkden(xx, x), lwd = 2, col = "red")
#' lines(density(x), lty = 2, lwd = 2, col = "green")
#' legend("topright", c("True Density", "KDE Using evmix", "KDE Using density function"),
#' lty = c(1, 1, 2), lwd = c(1, 2, 2), col = c("black", "red", "green"))
#' 
#' # Estimate bandwidth using cross-validation likelihood
#' x = rnorm(nk)
#' fit = fkden(x)
#' hist(x, nk/5, freq = FALSE, xlim = c(-5, 5), ylim = c(0, 0.6)) 
#' rug(x)
#' for (i in 1:nk) lines(xx, dnorm(xx, x[i], sd = fit$bw)*0.05)
#' lines(xx,dnorm(xx), col = "black")
#' lines(xx, dkden(xx, x, lambda = fit$lambda), lwd = 2, col = "red")
#' lines(density(x), lty = 2, lwd = 2, col = "green")
#' lines(density(x, bw = fit$bw), lwd = 2, lty = 2,  col = "blue")
#' legend("topright", c("True Density", "KDE fitted evmix",
#' "KDE Using density, default bandwidth", "KDE Using density, c-v likelihood bandwidth"),
#' lty = c(1, 1, 2, 2), lwd = c(1, 2, 2, 2), col = c("black", "red", "green", "blue"))
#'
#' plot(xx, pnorm(xx), type = "l")
#' rug(x)
#' lines(xx, pkden(xx, x), lwd = 2, col = "red")
#' lines(xx, pkden(xx, x, lambda = fit$lambda), lwd = 2, col = "green")
#' # green and blue (quantile) function should be same
#' p = seq(0, 1, 0.001)
#' lines(qkden(p, x, lambda = fit$lambda), p, lwd = 2, lty = 2, col = "blue") 
#' legend("topleft", c("True Density", "KDE using evmix, normal reference rule",
#' "KDE using evmix, c-v likelihood","KDE quantile function, c-v likelihood"),
#' lty = c(1, 1, 1, 2), lwd = c(1, 2, 2, 2), col = c("black", "red", "green", "blue"))
#' 
#' xnew = rkden(10000, x, lambda = fit$lambda)
#' hist(xnew, breaks = 100, freq = FALSE, xlim = c(-5, 5))
#' rug(xnew)
#' lines(xx,dnorm(xx), col = "black")
#' lines(xx, dkden(xx, x), lwd = 2, col = "red")
#' legend("topright", c("True Density", "KDE Using evmix"),
#' lty = c(1, 2), lwd = c(1, 2), col = c("black", "red"))
#' }
#' 
NULL

#' @export
#' @aliases kden dkden pkden qkden rkden
#' @rdname  kden

# density function for kernel density estimator
dkden <- function(x, kerncentres, lambda = NULL, bw = NULL, kernel = "gaussian", log = FALSE) {

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.quant(kerncentres, allowna = TRUE, allowinf = TRUE)
  check.kbw(lambda, bw, allownull = TRUE)
  check.kernel(kernel)
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

  d = x # pass through NA/NaN as entered
  
  whichok = which(is.finite(x))
  
  d[whichok] = sapply(x[whichok], FUN = kdenx, kerncentres = kerncentres,
    lambda = lambda, kernel = kernel) 

  if (log) d = log(d)

  d
}

#' @export
#' @aliases kden dkden pkden qkden rkden
#' @rdname  kden

# cumulative distribution function for kernel density estimator
pkden <- function(q, kerncentres, lambda = NULL, bw = NULL, kernel = "gaussian", lower.tail = TRUE) {

  # Check properties of inputs
  check.quant(q, allowna = TRUE, allowinf = TRUE)
  check.quant(kerncentres, allowna = TRUE, allowinf = TRUE)
  check.kbw(lambda, bw, allownull = TRUE)
  check.kernel(kernel)
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

  p = q # pass through NA/NaN as entered
  
  whichok = which(is.finite(q))
  
  p[whichok] = sapply(q[whichok], FUN = pkdenx, kerncentres = kerncentres,
    lambda = lambda, kernel = kernel) 

  # sometimes due to numerical errors p>1 or p<0
  p = pmax(pmin(p, 1), 0)
  
  if (!lower.tail) p = 1 - p

  p
}

#' @export
#' @aliases kden dkden pkden qkden rkden
#' @rdname  kden

# inverse cumulative distribution function for kernel density estimator
qkden <- function(p, kerncentres, lambda = NULL, bw = NULL, kernel = "gaussian", lower.tail = TRUE) {

  # Check properties of inputs
  check.prob(p, allowna = TRUE)
  check.quant(kerncentres, allowna = TRUE, allowinf = TRUE)
  check.kbw(lambda, bw, allownull = TRUE)
  check.kernel(kernel)
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

  if (!lower.tail) p = 1 - p

  q = p
  
  # obtain quantile function but interpolation between kernel based CDF estimates
  # - CDF is quick to estimate, interpolation algorithms also quite fast
  # - an alternative solution is to numerical solve to find quantile, but is much slower

  maxp = ifelse(kernel == "gaussian", 5*lambda, lambda)
  
  qrange = range(kerncentres) + maxp * c(-1, 1)
  qk = seq(qrange[1], qrange[2], length.out = min(50 * nk, 1000))

  pk = sapply(qk, FUN = pkdenx, lambda = lambda, kerncentres = kerncentres, kernel = kernel) 

  qfun = splinefun(x = pk, y = qk)

  whichinterp = which((p >= min(pk)) & (p <= max(pk)))
  ninterp = length(whichinterp)
  which0 = which(p < min(pk))
  n0 = length(which0)
  which1 = which(p > max(pk))
  n1 = length(which1)
  
  if (ninterp > 0) q[whichinterp] = qfun(p[whichinterp])

  # if further out than 5 standard deviations from kernel centres then set:
  # lower tail to -Inf and upper tail set to Inf
  if (n0 > 0) q[which0] = -Inf
  if (n1 > 0) q[which1] = Inf
  
  q
}

#' @export
#' @aliases kden dkden pkden qkden rkden
#' @rdname  kden

# random number generation for kernel density estimator
rkden <- function(n = 1, kerncentres, lambda = NULL, bw = NULL, kernel = "gaussian") {

  # Check properties of inputs
  check.n(n)
  check.quant(kerncentres, allowna = TRUE, allowinf = TRUE)
  check.kbw(lambda, bw, allownull = TRUE)
  check.kernel(kernel)

  # Inversion sampling (very inefficient!, but avoids need for RNG's for all different kernels)
  qkden(runif(n), kerncentres, lambda, bw, kernel)
}
