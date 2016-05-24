#' @name bckden
#' 
#' @title Boundary Corrected Kernel Density Estimation Using a Variety of Approaches
#'
#' @description Density, cumulative distribution function, quantile function and
#'   random number generation for boundary corrected kernel density estimators
#'   using a variety of approaches (and different kernels) with a constant
#'   bandwidth \code{lambda}.
#'
#' @param bcmethod boundary correction method
#' @param proper   logical, whether density is renormalised to integrate to unity (where needed)
#' @param nn       non-negativity correction method (simple boundary correction only)
#' @param xmax     upper bound on support (copula and beta kernels only) or \code{NULL}
#' @param offset   offset added to kernel centres (logtrans only) or \code{NULL}
#' @inheritParams kden
#' @inheritParams kernels
#' @inheritParams gpd
#'
#' @details Boundary corrected kernel density estimation (BCKDE) with improved
#' bias properties near the boundary compared to standard KDE available in 
#' \code{\link[evmix:kden]{kden}} functions. The user chooses from a wide range
#' of boundary correction methods designed to cope with a lower bound at zero
#' and potentially also both upper and lower bounds.
#' 
#' Some boundary correction methods require a secondary correction for
#' negative density estimates of which two methods are implemented. Further, some
#' methods don't necessarily give a density which integrates to one, so an option
#' is provided to renormalise to be proper.
#' 
#' It assumes there is a lower bound at zero, so prior transformation of data is
#' required for a alternative lower bound (possibly including negation to allow
#' for an upper bound).
#' 
#' The alternate bandwidth definitions are discussed in the
#' \code{\link[evmix:kernels]{kernels}}, with the \code{lambda} as the default.
#' The \code{bw} specification is the same as used in the
#' \code{\link[stats:density]{density}} function.
#' 
#' Certain boundary correction methods use the standard kernels which are defined
#' in the \code{\link[evmix:kernels]{kernels}} help
#' documentation with the \code{"gaussian"} as the default choice.
#'
#' The quantile function is rather complicated as there is no closed form solution,
#' so is obtained by numerical approximation of the inverse cumulative distribution function
#' \eqn{P(X \le q) = p} to find \eqn{q}. The quantile function 
#' \code{\link[evmix:bckden]{qbckden}} evaluates the KDE cumulative distribution
#' function over the range from \code{c(0, max(kerncentre) + lambda)},
#' or \code{c(0, max(kerncentre) + 5*lambda)} for normal kernel. Outside of this
#' range the quantiles are set to \code{0} for lower tail and \code{Inf}
#' (or \code{xmax} where appropriate) for upper tail. A sequence of values
#' of length fifty times the number of kernels (upto a maximum of 1000) is first
#' calculated. Spline based interpolation using \code{\link[stats:splinefun]{splinefun}},
#' with default \code{monoH.FC} method, is then used to approximate the quantile
#' function. This is a similar approach to that taken
#' by Matt Wand in the \code{\link[ks:kde.1d]{qkde}} in the \code{\link[ks:kde.1d]{ks}} package.
#' 
#' Unlike the standard KDE, there is no general rule-of-thumb bandwidth for all these
#' estimators, with only certain methods having a guideline in the literature, so none
#' have been implemented. Hence, a bandwidth must always be specified and you should
#' consider using \code{\link[evmix:fbckden]{fbckden}} function for cross-validation
#' MLE for bandwidth.
#' 
#' Random number generation is slow as inversion sampling using the (numerically evaluated)
#' quantile function is implemented. Users may want to consider alternative approaches instead,
#' like rejection sampling.
#' 
#' @section Boundary Correction Methods:
#' 
#' Renormalisation to a proper density is assumed by default \code{proper=TRUE}.
#' This correction is needed for \code{bcmethod="renorm"}, \code{"simple"},
#' \code{"beta1"}, \code{"beta2"}, \code{"gamma1"} and \code{"gamma2"} which
#' all require numerical integration. Renormalisation will not be carried out
#' for other methods, even when \code{proper=TRUE}.
#' 
#' Non-negativity correction is only relevant for the \code{bcmethod="simple"} approach.
#' The Jones and Foster (1996) method is applied \code{nn="jf96"} by default. This method
#' can occassionally give an extra boundary bias for certain populations (e.g. Gamma(2, 1)),
#' see paper for details. Non-negative values can simply be zeroed (\code{nn="zero"}).
#' Renormalisation should always be applied after non-negativity correction. Non-negativity
#' correction will not be carried out for other methods, even when requested by user.
#' 
#' The non-negative correction is applied before renormalisation, when both requested. 
#' 
#' The boundary correction methods implemented are listed below. The first set can use
#' any type of kernel (see \code{\link[evmix:kernels]{kernels}} help
#' documentation):
#' 
#' \code{bcmethod="simple"} is the default and applies the simple boundary correction method
#' in equation (3.4) of Jones (1993) and is equivalent to the kernel weighted local linear
#' fitting at the boundary. Renormalisation and non-negativity correction may be required.
#' 
#' \code{bcmethod="cutnorm"} applies cut and normalisation method of
#' Gasser and Muller (1979), where the kernels themselves are individually truncated at
#' the boundary and renormalised to unity.
#' 
#' \code{bcmethod="renorm"} applies first order correction method discussed in
#' Diggle (1985), where the kernel density estimate is locally renormalised near boundary.
#' Renormalisation may be required.
#' 
#' \code{bcmethod="reflect"} applies reflection method of Boneva, Kendall and Stefanov
#' (1971) which is equivalent to the dataset being supplemented by the same dataset negated. 
#' This method implicitly assumes f'(0)=0, so can cause extra artefacts at the boundary. 
#' 
#' \code{bcmethod="logtrans"} applies KDE on the log-scale and then back-transforms (with
#' explicit normalisation) following Marron and Ruppert (1992). This is the approach
#' implemented in the \code{\link[ks:kde.1d]{ks}} package. As the KDE is applied on
#' the log scale, the effective bandwidth on the original scale is non-constant. The
#' \code{offset} option is only used for this method and is commonly used to offset
#' zero kernel centres in log transform to prevent \code{log(0)}.
#' 
#' All the following boundary correction methods do not use kernels in their
#' usual sense, so ignore the \code{kernel} input:
#' 
#' \code{bcmethod="beta1"} and \code{"beta2"} uses the beta and modified beta kernels
#' of Chen (1999) respectively. The \code{xmax} rescales the beta kernels to be
#' defined on the support [0, xmax] rather than unscaled [0, 1]. Renormalisation
#' will be required.
#' 
#' \code{bcmethod="gamma1"} and \code{"gamma2"} uses the gamma and modified gamma kernels
#' of Chen (2000) respectively. Renormalisation will be required.
#' 
#' \code{bcmethod="copula"} uses the bivariate normal copula based kernesl of 
#' Jones and Henderson (2007). As with the \code{bcmethod="beta1"}  and \code{"beta2"}
#' methods the \code{xmax} rescales the copula kernels to be defined on the support [0, xmax]
#' rather than [0, 1]. In this case the bandwidth is defined as \eqn{lambda=1-\rho^2},
#' so the bandwidth is limited to \eqn{(0, 1)}.
#' 
#' @section Warning:
#' The \code{"simple"}, \code{"renorm"}, \code{"beta1"}, \code{"beta2"}, \code{"gamma1"} 
#' and \code{"gamma2"} boundary correction methods may require renormalisation using
#' numerical integration which can be very slow. In particular, the numerical integration
#' is extremely slow for the \code{kernel="uniform"}, due to the adaptive quadrature in
#' the \code{\link[stats:integrate]{integrate}} function
#' being particularly slow for functions with step-like behaviour.
#' 
#' @return \code{\link[evmix:bckden]{dbckden}} gives the density, 
#' \code{\link[evmix:bckden]{pbckden}} gives the cumulative distribution function,
#' \code{\link[evmix:bckden]{qbckden}} gives the quantile function and 
#' \code{\link[evmix:bckden]{rbckden}} gives a random sample.
#' 
#' @note Unlike most of the other extreme value mixture model functions the 
#' \code{\link[evmix:bckden]{bckden}} functions have not been vectorised as
#' this is not appropriate. The main inputs (\code{x}, \code{p} or \code{q})
#' must be either a scalar or a vector, which also define the output length.
#' 
#' The kernel centres \code{kerncentres} can either be a single datapoint or a vector
#' of data. The kernel centres (\code{kerncentres}) and locations to evaluate density (\code{x})
#' and cumulative distribution function (\code{q}) would usually be different.
#' 
#' Default values are provided for all inputs, except for the fundamentals 
#' \code{lambda}, \code{kerncentres}, \code{x}, \code{q} and \code{p}.
#' The default sample size for \code{\link[evmix:bckden]{rbckden}} is 1.
#' 
#' The \code{xmax} option is only relevant for the beta and copula methods, so a
#' warning is produced if this is not \code{NULL} for in other methods.
#' The \code{offset} option is only relevant for the \code{"logtrans"} method, so a
#' warning is produced if this is not \code{NULL} for in other methods.
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
#' Chen, S.X. (1999). Beta kernel estimators for density functions. Computational Statistics
#' and Data Analysis 31, 1310-45.
#' 
#' Gasser, T. and Muller, H. (1979). Kernel estimation of regression functions. In "Lecture Notes
#' in Mathematics 757, edited by Gasser and Rosenblatt, Springer.
#' 
#' Chen, S.X. (2000). Probability density function estimation using gamma kernels.
#' Annals of the Institute of Statisical Mathematics 52(3), 471-480.
#' 
#' Boneva, L.I., Kendall, D.G. and Stefanov, I. (1971). Spline transformations: Three new
#' diagnostic aids for the statistical data analyst (with discussion). Journal of the Royal
#' Statistical Society B, 33, 1-70.
#' 
#' Diggle, P.J. (1985). A kernel method for smoothing point process data. Applied Statistics
#' 34, 138-147.
#' 
#' Marron, J.S. and Ruppert, D. (1994) Transformations to reduce boundary bias in kernel
#' density estimation, Journal of the Royal Statistical Society. Series B 56(4), 653-671.
#' 
#' Jones, M.C. and Henderson, D.A. (2007). Kernel-type density estimation on the unit
#' interval. Biometrika 94(4), 977-984.
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
#' @aliases bckden dbckden pbckden qbckden rbckden
#' @family  kden kdengpd kdengpdcon bckden bckdengpd bckdengpdcon
#'          fkden fkdengpd fkdengpdcon fbckden fbckdengpd fbckdengpdcon
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' par(mfrow = c(1, 1))
#' 
#' n=100
#' x = rgamma(n, shape = 1, scale = 2)
#' xx = seq(-0.5, 12, 0.01)
#' plot(xx, dgamma(xx, shape = 1, scale = 2), type = "l")
#' rug(x)
#' lines(xx, dbckden(xx, x, lambda = 1), lwd = 2, col = "red")
#' lines(density(x), lty = 2, lwd = 2, col = "green")
#' legend("topright", c("True Density", "Simple boundary correction",
#' "KDE using density function", "Boundary Corrected Kernels"),
#' lty = c(1, 1, 2, 1), lwd = c(1, 2, 2, 1), col = c("black", "red", "green", "blue"))
#'
#' n=100
#' x = rbeta(n, shape1 = 3, shape2 = 2)*5
#' xx = seq(-0.5, 5.5, 0.01)
#' plot(xx, dbeta(xx/5, shape1 = 3, shape2 = 2)/5, type = "l", ylim = c(0, 0.8))
#' rug(x)
#' lines(xx, dbckden(xx, x, lambda = 0.1, bcmethod = "beta2", proper = TRUE, xmax = 5),
#'   lwd = 2, col = "red")
#' lines(density(x), lty = 2, lwd = 2, col = "green")
#' legend("topright", c("True Density", "Modified Beta KDE Using evmix",
#'   "KDE using density function"),
#' lty = c(1, 1, 2), lwd = c(1, 2, 2), col = c("black", "red", "green"))
#'
#' # Demonstrate renormalisation (usually small difference)
#' n=1000
#' x = rgamma(n, shape = 1, scale = 2)
#' xx = seq(-0.5, 15, 0.01)
#' plot(xx, dgamma(xx, shape = 1, scale = 2), type = "l")
#' rug(x)
#' lines(xx, dbckden(xx, x, lambda = 0.5, bcmethod = "simple", proper = TRUE),
#'   lwd = 2, col = "purple")
#' lines(xx, dbckden(xx, x, lambda = 0.5, bcmethod = "simple", proper = FALSE),
#'   lwd = 2, col = "red", lty = 2)
#' legend("topright", c("True Density", "Simple BC with renomalisation", 
#' "Simple BC without renomalisation"),
#' lty = 1, lwd = c(1, 2, 2), col = c("black", "purple", "red"))
#' }
#' 
NULL

#' @export
#' @aliases bckden dbckden pbckden qbckden rbckden
#' @rdname  bckden

# density function for boundary corrected KDE
dbckden <- function(x, kerncentres, lambda = NULL, bw = NULL, kernel = "gaussian",
  bcmethod = "simple", proper = TRUE, nn = "jf96", offset = NULL, xmax = NULL, log = FALSE) {

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.quant(kerncentres, allowna = TRUE, allowinf = TRUE)
  check.kbw(lambda, bw)
  check.kernel(kernel)
  check.bcmethod(bcmethod)
  check.logic(proper)
  check.nn(nn)
  check.offset(offset, bcmethod, allowzero = TRUE)
  check.posparam(xmax, allownull = TRUE)  
  check.logic(log)

  kernel = ifelse(kernel == "rectangular", "uniform", kernel)
  kernel = ifelse(kernel == "normal", "gaussian", kernel)

  if (any(is.infinite(x))) warning("infinite quantiles set to NA")

  x[is.infinite(x)] = NA # user will have to deal with infinite cases
    
  if (any(!is.finite(kerncentres))) warning("non-finite kernel centres are dropped")

  kerncentres = kerncentres[is.finite(kerncentres)]

  if (any(kerncentres < 0)) stop("negative kernel centres not permitted")

  check.quant(kerncentres)
  nk = length(kerncentres)
  
  # if bcmethod does not use standard kernels then lambda must be specified
  # then bw can be used, but lambda should be defaulted to if available
  kernelmethods = c("simple", "cutnorm", "renorm", "reflect", "logtrans")
  if (!(bcmethod %in% kernelmethods)) {
    if (is.null(lambda))
      stop(paste("bandwidth bw only relevant for", kernelmethods, collapse = " "))
  } else {
    lambda = klambda(bw, kernel, lambda)    
  }

  if ((bcmethod == "copula") & (lambda >= 1))
    stop("bandwidth must between (0, 1) for copula method")  
  
  upboundmethods = c("beta1", "beta2", "copula")
  if (!is.null(xmax) & !(bcmethod %in% upboundmethods))
    warning(paste("xmax only relevant for boundary correction methods", upboundmethods, collapse = " "))
  
  if (bcmethod %in% upboundmethods) {
    if (is.null(xmax)) stop("xmax is NULL")
    
    if (max(kerncentres) > xmax) stop("largest kernel centre must be below xmax")

    if (any(kerncentres == 0)) {
      warning("kernel centres of zero are invalid for gamma or beta method so ignored")
      kerncentres = kerncentres[kerncentres != 0]
    }

    if ((bcmethod != "gamma1") & (bcmethod != "gamma2")) {
      if (any(kerncentres == xmax)) {
        warning("kernel centres of xmax are invalid for beta or copula method so ignored")
        kerncentres = kerncentres[kerncentres != xmax]
      }
    }
    # need to recheck there are some valid kernel centres after these exclusions
    check.quant(kerncentres)
  }
  
  # numerical integration can be problematic if no data near boundary
  # bounds the evaluation range (zero outside of [minaccept, maxaccept])
  maxp = ifelse(kernel == "gaussian", 5, 1)*lambda
  if (bcmethod %in% upboundmethods) {
    minaccept = 0
    maxaccept = xmax
  } else if (bcmethod == "logtrans") {
    maxaccept = exp(log(max(kerncentres) + offset) + maxp)
    minaccept = max(offset + .Machine$double.eps*2, exp(log(min(kerncentres) + offset) - maxp))
  } else if (bcmethod == "gamma1") {
    maxaccept = qgamma(1e-8, shape = max(kerncentres)/lambda + 1, scale = lambda, lower.tail = FALSE)
    minaccept = qgamma(1e-8, shape = min(kerncentres)/lambda + 1, scale = lambda)
  } else if (bcmethod == "gamma2") {
    maxaccept = ifelse(max(kerncentres) > 2*lambda, 
      qgamma(1e-8, shape = max(kerncentres)/lambda, scale = lambda, lower.tail = FALSE),
      qgamma(1e-8, shape = (max(kerncentres)/lambda)^2/4 + 1, scale = lambda, lower.tail = FALSE))
    minaccept = ifelse(min(kerncentres) > 2*lambda, 
      qgamma(1e-8, shape = min(kerncentres)/lambda, scale = lambda),
      qgamma(1e-8, shape = (min(kerncentres)/lambda)^2/4 + 1, scale = lambda))
  } else {
    maxaccept = max(kerncentres) + maxp
    minaccept = max(0, min(kerncentres) - maxp)
  }
  
  d = x # pass through NA/NaN as entered
  
  d[!is.na(d)] = 0 # default outside of support

  # only select non-negative non-missing x-values for evaluation
  xok = x[!is.na(x)]
  xok = xok[(xok >= minaccept) & (xok <= maxaccept)]
  nok = length(xok)
  
  if (nok > 0) {
    if (bcmethod == "simple") {
      # simple linear boundary correction method of Jones (1993), eq 3.4
      # which adapts kernels to given (kernel weighted) local linear fit near boundary
      
      dok = sapply(xok, FUN = bckdenxsimple, kerncentres = kerncentres, lambda = lambda, kernel = kernel)
      
    } else if (bcmethod == "cutnorm") {
      # cut (at boundary) and normalisation each individual kernel
      
      dok = sapply(xok, FUN = bckdenxcutnorm, kerncentres = kerncentres, lambda = lambda, kernel = kernel)
      proper = FALSE # not needed
      nn = "none"    # not needed
      
    } else if (bcmethod == "renorm") {
      # first order correction of KDE near boundary
      
      dok = sapply(xok, FUN = bckdenxrenorm, kerncentres = kerncentres, lambda = lambda, kernel = kernel)
      nn = "none"    # not needed
      
    } else if (bcmethod == "reflect") {
      # reflection method which is equivalent to using (kerncentres, -kerncentres)
      # only good if f'(0)=0
      
      dok = sapply(xok, FUN = bckdenxreflect, kerncentres = kerncentres, lambda = lambda, kernel = kernel)
      proper = FALSE # not needed
      nn = "none"    # not needed
      
    } else if (bcmethod == "logtrans") {
      # log transformation with standard KDE applied on log-scale
      # (Wand, Marron and Ruppert JASA, 1991, Marron and Ruppert JRSS B, 1994)
      # requires offset in min(kerncentres) = 0
      
      # transformation KDE is f(x) = g(t(x)) abs(t'(x))
      dok = sapply(log(xok + offset), FUN = kdenx, 
        kerncentres = log(kerncentres + offset), lambda = lambda, kernel = kernel) / (xok + offset)
      proper = FALSE # not needed
      nn = "none"    # not needed
      
    } else if (bcmethod == "beta1") {
      # beta kernels by Chen (1999) CSDA
      
      dok = sapply(xok, FUN = bckdenxbeta1, kerncentres = kerncentres, lambda = lambda, xmax = xmax)
      nn = "none"    # not needed
      
    } else if (bcmethod == "beta2") {
      # modified beta kernels by Chen (1999) CSDA
      
      dok = sapply(xok, FUN = bckdenxbeta2, kerncentres = kerncentres, lambda = lambda, xmax = xmax)
      nn = "none"    # not needed
      
    } else if (bcmethod == "copula") {
      # bivariate normal copula based kernels by Jones and Henderson (2007) Biometrika
      
      dok = sapply(xok, FUN = bckdenxcopula, kerncentres = kerncentres, lambda = lambda, xmax = xmax)
      proper = FALSE # not needed
      nn = "none"    # not needed
      
    } else if (bcmethod == "gamma1") {
      # gamma kernels by Chen (2000) AISM
        
      dok = sapply(xok, FUN = bckdenxgamma1, kerncentres = kerncentres, lambda = lambda)
      nn = "none"    # not needed
      
    } else if (bcmethod == "gamma2") {
      # modified gamma kernels by Chen (2000) AISM
      
      dok = sapply(xok, FUN = bckdenxgamma2, kerncentres = kerncentres, lambda = lambda)
      nn = "none"    # not needed
      
    }
    
    # Apply non-negative correction first followed by renormalisation (if either requested)
    if (nn == "jf96") {
      dbar = sapply(xok, FUN = bckdenxrenorm, kerncentres = kerncentres, lambda = lambda, kernel = kernel)
      dok = ifelse(dbar == 0, dok, dbar*exp(dok/dbar - 1))
    } else if (nn =="zero") {
      dok[which(dok < 0)] = 0
    }

    if (proper) {
      if ((bcmethod == "simple") & (nn == "none")) {
        pxmax = pbckdenxsimple(maxaccept, kerncentres = kerncentres, lambda = lambda, kernel = kernel)
      } else if ((bcmethod == "simple") & (nn != "none")) {
        pxmax = pbckdenxnn(maxaccept, kerncentres = kerncentres, lambda = lambda, kernel = kernel, nn = nn)
      } else if (bcmethod == "renorm") {
        pxmax = pbckdenxrenorm(maxaccept, kerncentres = kerncentres, lambda = lambda, kernel = kernel)
      } else if (bcmethod == "beta1") {
        pxmax = pbckdenxbeta1(maxaccept, kerncentres = kerncentres, lambda = lambda, xmax = xmax)
      } else if (bcmethod == "beta2") {
        pxmax = pbckdenxbeta2(maxaccept, kerncentres = kerncentres, lambda = lambda, xmax = xmax)
      } else if (bcmethod == "gamma1") {
        pxmax = pbckdenxgamma1(maxaccept, kerncentres = kerncentres, lambda = lambda)
      } else if (bcmethod == "gamma2") {
        pxmax = pbckdenxgamma2(maxaccept, kerncentres = kerncentres, lambda = lambda)
      }
      dok = dok/pxmax
    }
    d[ifelse(!is.na(x), (x >= minaccept) & (x <= maxaccept), FALSE)] = dok
  }
  
  if (log) d = log(d)

  d
}


#' @export
#' @aliases bckden dbckden pbckden qbckden rbckden
#' @rdname  bckden

# cumulative distribution function for boundary corrected KDE
pbckden <- function(q, kerncentres, lambda = NULL, bw = NULL, kernel = "gaussian",
  bcmethod = "simple", proper = TRUE, nn = "jf96", offset = NULL, xmax = NULL, lower.tail = TRUE) {

  # Check properties of inputs
  check.quant(q, allowna = TRUE, allowinf = TRUE)
  check.quant(kerncentres, allowna = TRUE, allowinf = TRUE)
  check.kbw(lambda, bw)
  check.kernel(kernel)
  check.bcmethod(bcmethod)
  check.logic(proper)
  check.nn(nn)
  check.offset(offset, bcmethod, allowzero = TRUE)
  check.posparam(xmax, allownull = TRUE)  
  check.logic(lower.tail)

  kernel = ifelse(kernel == "rectangular", "uniform", kernel)
  kernel = ifelse(kernel == "normal", "gaussian", kernel)

  if (any(is.infinite(q))) warning("infinite quantiles set to NA")

  q[is.infinite(q)] = NA # user will have to deal with infinite cases
    
  if (any(!is.finite(kerncentres))) warning("non-finite kernel centres are dropped")

  kerncentres = kerncentres[is.finite(kerncentres)]

  if (any(kerncentres < 0)) stop("negative kernel centres not permitted")

  check.quant(kerncentres)
  nk = length(kerncentres)
  
  # if bcmethod does not use standard kernels then lambda must be specified
  # then bw can be used, but lambda should be defaulted to if available
  kernelmethods = c("simple", "cutnorm", "renorm", "reflect", "logtrans")
  if (!(bcmethod %in% kernelmethods)) {
    if (is.null(lambda))
      stop(paste("bandwidth bw only relevant for", kernelmethods, collapse = " "))
  } else {
    lambda = klambda(bw, kernel, lambda)    
  }

  if ((bcmethod == "copula") & (lambda >= 1))
    stop("bandwidth must between (0, 1) for copula method")  
    
  upboundmethods = c("beta1", "beta2", "copula")
  if (!is.null(xmax) & !(bcmethod %in% upboundmethods))
    warning(paste("xmax only relevant for boundary correction methods", upboundmethods, collapse = " "))
  
  if (bcmethod %in% upboundmethods) {
    if (is.null(xmax)) stop("xmax is NULL")
    
    if (max(kerncentres) > xmax) stop("largest kernel centre must be below xmax")

    if (any(kerncentres == 0)) {
      warning("kernel centres of zero are invalid for gamma or beta method so ignored")
      kerncentres = kerncentres[kerncentres != 0]
    }

    if ((bcmethod != "gamma1") & (bcmethod != "gamma2")) {
      if (any(kerncentres == xmax)) {
        warning("kernel centres of xmax are invalid for beta or copula method so ignored")
        kerncentres = kerncentres[kerncentres != xmax]
      }
    }
    # need to recheck there are some valid kernel centres after these exclusions
    check.quant(kerncentres)
  }
   
  # numerical integration can be problematic if no data near boundary
  # bounds the evaluation range (zero outside of [minaccept, maxaccept])
  maxp = ifelse(kernel == "gaussian", 5, 1)*lambda
  if (bcmethod %in% upboundmethods) {
    minaccept = 0
    maxaccept = xmax
  } else if (bcmethod == "logtrans") {
    maxaccept = exp(log(max(kerncentres) + offset) + maxp)
    minaccept = max(offset + .Machine$double.eps*2, exp(log(min(kerncentres) + offset) - maxp))
  } else if (bcmethod == "gamma1") {
    maxaccept = qgamma(1e-6, shape = max(kerncentres)/lambda + 1, scale = lambda, lower.tail = FALSE)
    minaccept = qgamma(1e-6, shape = min(kerncentres)/lambda + 1, scale = lambda)
  } else if (bcmethod == "gamma2") {
    maxaccept = ifelse(max(kerncentres) > 2*lambda, 
      qgamma(1e-6, shape = max(kerncentres)/lambda, scale = lambda, lower.tail = FALSE),
      qgamma(1e-6, shape = (max(kerncentres)/lambda)^2/4 + 1, scale = lambda, lower.tail = FALSE))
    minaccept = ifelse(min(kerncentres) > 2*lambda, 
      qgamma(1e-6, shape = min(kerncentres)/lambda, scale = lambda),
      qgamma(1e-6, shape = (min(kerncentres)/lambda)^2/4 + 1, scale = lambda))
  } else {
    maxaccept = max(kerncentres) + maxp
    minaccept = max(0, min(kerncentres) - maxp)
  }

  p = q # pass through NA/NaN as entered

  p[!is.na(q)] = ifelse(q[!is.na(q)] < minaccept, 0, p)

  # if normalised then maximum value is 1
  # if not normalised then integrate upto maxaccept to find upper limit
  if (proper) {
    p[which(!is.na(q))] = ifelse(q[which(!is.na(q))] > maxaccept, 1, p)
  } else {
    p[which(!is.na(q))] = ifelse(q[which(!is.na(q))] > maxaccept, 
      pbckden(maxaccept, kerncentres, lambda, kernel = kernel, 
        bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax), p)
  }
  
  # only select non-negative x-values for evaluation
  qok = q[!is.na(q)]
  qok = qok[(qok >= minaccept) & (qok <= maxaccept)]
  nok = length(qok)
  
  if (nok > 0) {
    if ((bcmethod == "simple") & (nn == "none")) {
      # simple linear boundary correction method of Jones (1993), eq 3.4
      # which adapts kernels to given (kernel weighted) local linear fit near boundary
      # no negative-boundary correction
      
      pok = sapply(qok, FUN = pbckdenxsimple, kerncentres = kerncentres, lambda = lambda, kernel = kernel)
            
    } else if ((bcmethod == "simple") & (nn != "none")) {
      # same as above, but also applies non-negative density correction
      # of Jones and Foster (1996) or crude zeroing
      
      pok = sapply(qok, FUN = pbckdenxnn, kerncentres = kerncentres, lambda = lambda, kernel = kernel, nn = nn)
      
    } else if (bcmethod == "cutnorm") {
      # cut (at boundary) and normalisation each individual kernel
      
      pok = sapply(qok, FUN = pbckdenxcutnorm, kerncentres = kerncentres, lambda = lambda, kernel = kernel)
      pok = ifelse(pok > 1, 1, pok) # numerical errors occassionally give cdf > 1
      proper = FALSE # not needed
      nn = "none"    # not needed
      
    } else if (bcmethod == "renorm") {
      # first order correction of KDE near boundary
      
      pok = sapply(qok, FUN = pbckdenxrenorm, kerncentres = kerncentres, lambda = lambda, kernel = kernel)
      nn = "none"    # not needed
      
    } else if (bcmethod == "reflect") {
      # reflection method which is equivalent to using (kerncentres, -kerncentres)
      # only good if f'(0)=0
            
      pok = sapply(qok, FUN = pbckdenxreflect, kerncentres = kerncentres, lambda = lambda, kernel = kernel)
      pok = ifelse(pok > 1, 1, pok) # numerical errors occassionally give cdf > 1
      proper = FALSE # not needed
      nn = "none"    # not needed
      
    } else if (bcmethod == "logtrans") {
      # log transformation with standard KDE applied on log-scale
      # (Wand, Marron and Ruppert JASA, 1991, Marron and Ruppert JRSS B, 1994)
      # requires offset in min(kerncentres) = 0
      
      # transformation KDE is f(x) = g(t(x)) abs(t'(x))
      pok = sapply(qok, FUN = pbckdenxlog, kerncentres = kerncentres, lambda = lambda,
        offset = offset, kernel = kernel)
      pok = ifelse(pok > 1, 1, pok) # numerical errors occassionally give cdf > 1
      proper = FALSE # not needed
      nn = "none"    # not needed
      
    } else if (bcmethod == "beta1") {
      # beta kernels by Chen (1999) CSDA
      
      pok = sapply(qok, FUN = pbckdenxbeta1, kerncentres = kerncentres, lambda = lambda, xmax = xmax)
      nn = "none"    # not needed
      
    } else if (bcmethod == "beta2") {
      # modified beta kernels by Chen (1999) CSDA
      
      pok = sapply(qok, FUN = pbckdenxbeta2, kerncentres = kerncentres, lambda = lambda, xmax = xmax)
      nn = "none"    # not needed
      
    } else if (bcmethod == "copula") {
      # bivariate normal copula based kernels by Jones and Henderson (2007) Biometrika
      
      pok = sapply(qok, FUN = pbckdenxcopula, kerncentres = kerncentres, lambda = lambda, xmax = xmax)
      pok = ifelse(pok > 1, 1, pok) # numerical errors occassionally give cdf > 1
      proper = FALSE # not needed
      nn = "none"    # not needed

    } else if (bcmethod == "gamma1") {
      # gamma kernels by Chen (2000) AISM

      pok = sapply(qok, FUN = pbckdenxgamma1, kerncentres = kerncentres, lambda = lambda)
      nn = "none"    # not needed
     
    } else if (bcmethod == "gamma2") {
      # modified gamma kernels by Chen (2000) AISM
      
      pok = sapply(qok, FUN = pbckdenxgamma2, kerncentres = kerncentres, lambda = lambda)
      nn = "none"    # not needed
      
    }
    
    if (proper) {
      if ((bcmethod == "simple") & (nn == "none")) {
        pxmax = pbckdenxsimple(maxaccept, kerncentres = kerncentres, lambda = lambda, kernel = kernel)
      } else if ((bcmethod == "simple") & (nn != "none")) {
        pxmax = pbckdenxnn(maxaccept, kerncentres = kerncentres, lambda = lambda, kernel = kernel, nn = nn)
      } else if (bcmethod == "renorm") {
        pxmax = pbckdenxrenorm(maxaccept, kerncentres = kerncentres, lambda = lambda, kernel = kernel)
      } else if (bcmethod == "beta1") {
        pxmax = pbckdenxbeta1(maxaccept, kerncentres = kerncentres, lambda = lambda, xmax = xmax)
      } else if (bcmethod == "beta2") {
        pxmax = pbckdenxbeta2(maxaccept, kerncentres = kerncentres, lambda = lambda, xmax = xmax)
      } else if (bcmethod == "gamma1") {
        pxmax = pbckdenxgamma1(maxaccept, kerncentres = kerncentres, lambda = lambda)
      } else if (bcmethod == "gamma2") {
        pxmax = pbckdenxgamma2(maxaccept, kerncentres = kerncentres, lambda = lambda)
      }
      pok = pok/pxmax
      pok = ifelse(pok > 1, 1, pok) # numerical errors occassionally give cdf > 1
    }
    
    p[ifelse(!is.na(q), (q >= minaccept) & (q <= maxaccept), FALSE)] = pok
  }
  
  if (!lower.tail) p = 1 - p
  p
}

#' @export
#' @aliases bckden dbckden pbckden qbckden rbckden
#' @rdname  bckden

# inverse cumulative distribution function for boundary corrected KDE
qbckden <- function(p, kerncentres, lambda = NULL, bw = NULL, kernel = "gaussian",
  bcmethod = "simple", proper = TRUE, nn = "jf96", offset = NULL, xmax = NULL, lower.tail = TRUE) {

  # Check properties of inputs
  check.prob(p, allowna = TRUE)
  check.quant(kerncentres, allowna = TRUE, allowinf = TRUE)
  check.kbw(lambda, bw)
  check.kernel(kernel)
  check.bcmethod(bcmethod)
  check.logic(proper)
  check.nn(nn)
  check.offset(offset, bcmethod, allowzero = TRUE)
  check.posparam(xmax, allownull = TRUE)
  check.logic(lower.tail)
    
  kernel = ifelse(kernel == "rectangular", "uniform", kernel)
  kernel = ifelse(kernel == "normal", "gaussian", kernel)

  if (any(!is.finite(kerncentres))) warning("non-finite kernel centres are dropped")

  kerncentres = kerncentres[is.finite(kerncentres)]

  if (any(kerncentres < 0)) stop("negative kernel centres not permitted")

  check.quant(kerncentres)
  nk = length(kerncentres)
  
  # if bcmethod does not use standard kernels then lambda must be specified
  # then bw can be used, but lambda should be defaulted to if available
  kernelmethods = c("simple", "cutnorm", "renorm", "reflect", "logtrans")
  if (!(bcmethod %in% kernelmethods)) {
    if (is.null(lambda))
      stop(paste("bandwidth bw only relevant for", kernelmethods, collapse = " "))
  } else {
    lambda = klambda(bw, kernel, lambda)    
  }

  if ((bcmethod == "copula") & (lambda >= 1))
    stop("bandwidth must between (0, 1) for copula method")  
    
  upboundmethods = c("beta1", "beta2", "copula")
  if (!is.null(xmax) & !(bcmethod %in% upboundmethods))
    warning(paste("xmax only relevant for boundary correction methods", upboundmethods, collapse = " "))
  
  if (bcmethod %in% upboundmethods) {
    if (is.null(xmax)) stop("xmax is NULL")
    
    if (max(kerncentres) > xmax) stop("largest kernel centre must be below xmax")

    if (any(kerncentres == 0)) {
      warning("kernel centres of zero are invalid for gamma or beta method so ignored")
      kerncentres = kerncentres[kerncentres != 0]
    }

    if ((bcmethod != "gamma1") & (bcmethod != "gamma2")) {
      if (any(kerncentres == xmax)) {
        warning("kernel centres of xmax are invalid for beta or copula method so ignored")
        kerncentres = kerncentres[kerncentres != xmax]
      }
    }
    # need to recheck there are some valid kernel centres after these exclusions
    check.quant(kerncentres)
  }
  
  # numerical integration can be problematic if no data near boundary
  # bounds the evaluation range (zero outside of [minaccept, maxaccept])
  maxp = ifelse(kernel == "gaussian", 5, 1)*lambda
  if (bcmethod %in% upboundmethods) {
    minaccept = 0
    maxaccept = xmax
  } else if (bcmethod == "logtrans") {
    maxaccept = exp(log(max(kerncentres) + offset) + maxp)
    minaccept = max(offset + .Machine$double.eps*2, exp(log(min(kerncentres) + offset) - maxp))
  } else if (bcmethod == "gamma1") {
    maxaccept = qgamma(1e-6, shape = max(kerncentres)/lambda + 1, scale = lambda, lower.tail = FALSE)
    minaccept = qgamma(1e-6, shape = min(kerncentres)/lambda + 1, scale = lambda)
  } else if (bcmethod == "gamma2") {
    maxaccept = ifelse(max(kerncentres) > 2*lambda, 
      qgamma(1e-6, shape = max(kerncentres)/lambda, scale = lambda, lower.tail = FALSE),
      qgamma(1e-6, shape = (max(kerncentres)/lambda)^2/4 + 1, scale = lambda, lower.tail = FALSE))
    minaccept = ifelse(min(kerncentres) > 2*lambda, 
      qgamma(1e-6, shape = min(kerncentres)/lambda, scale = lambda),
      qgamma(1e-6, shape = (min(kerncentres)/lambda)^2/4 + 1, scale = lambda))
  } else {
    maxaccept = max(kerncentres) + maxp
    minaccept = max(0, min(kerncentres) - maxp)
  }

  if (!lower.tail) p = 1 - p
  
  q = p # pass through NA/NaN as entered
  
  pok = p[!is.na(p)]
  
  qok = rep(NA, length(pok)) # initialise
  
  qk = seq(minaccept, maxaccept, length.out = ifelse(bcmethod == "logtrans", 1000, min(50*nk, 1000)))

  pk = sapply(qk, FUN = pbckden, kerncentres = kerncentres, lambda = lambda, kernel = kernel,
    bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax)

  qfun = splinefun(x = pk, y = qk, method = "monoH.FC") # good method for cdf as monotone

  whichinterp = which((pok > min(pk)) & (pok < max(pk)))
  ninterp = length(whichinterp)
  which0 = which(pok <= min(pk))
  n0 = length(which0)
  which1 = which(pok >= max(pk))
  n1 = length(which1)

  if (ninterp > 0) qok[whichinterp] = qfun(pok[whichinterp])

  # if further out then set:
  # lower tail to -Inf and upper tail set to Inf
  if (n0 > 0) qok[which0] = 0
  if (n1 > 0) qok[which1] = rep(ifelse(bcmethod %in% upboundmethods, xmax, Inf), length(which1))
  
  q[!is.na(p)] = qok

  q
}

#' @export
#' @aliases bckden dbckden pbckden qbckden
#' @rdname  bckden

# random number generation for boundary corrected KDE
rbckden <- function(n = 1, kerncentres, lambda = NULL, bw = NULL, kernel = "gaussian",
  bcmethod = "simple", proper = TRUE, nn = "jf96", offset = NULL, xmax = NULL) {

  # Check properties of inputs
  check.n(n)
  check.quant(kerncentres, allowna = TRUE, allowinf = TRUE)
  check.kbw(lambda, bw)
  check.kernel(kernel)
  check.bcmethod(bcmethod)
  check.logic(proper)
  check.nn(nn)
  check.offset(offset, bcmethod, allowzero = TRUE)
  check.posparam(xmax, allownull = TRUE)  

  qbckden(runif(n), kerncentres, lambda, bw, kernel, bcmethod, proper, nn, offset, xmax)
}
