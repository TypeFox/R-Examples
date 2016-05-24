#' @name bckdengpd
#' 
#' @title Boundary Corrected Kernel Density Estimate and GPD Tail Extreme Value Mixture Model
#'
#' @description Density, cumulative distribution function, quantile function and
#'   random number generation for the extreme value mixture model with 
#'   boundary corrected kernel density estimate for bulk
#'   distribution upto the threshold and conditional GPD above threshold. The parameters
#'   are the bandwidth \code{lambda}, threshold \code{u}
#'   GPD scale \code{sigmau} and shape \code{xi} and tail fraction \code{phiu}.
#'
#' @inheritParams bckden
#' @inheritParams kdengpd
#' @inheritParams gpd
#' 
#' @details Extreme value mixture model combining boundary corrected kernel density (BCKDE)
#' estimate for the bulk below the threshold and GPD for upper tail. The user chooses
#' from a wide range of boundary correction methods designed to cope with a lower bound
#' at zero and potentially also both upper and lower bounds.
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
#' The user can pre-specify \code{phiu} permitting a parameterised value for the
#' tail fraction \eqn{\phi_u}. Alternatively, when \code{phiu=TRUE} the tail fraction
#' is estimated as the tail fraction from the BCKDE bulk model.
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
#' upper tail fraction of the BCKDE (\code{phiu=TRUE}), upto the threshold
#' \eqn{x \le u}, given by:
#' \deqn{F(x) = H(x)}
#' and above the threshold \eqn{x > u}:
#' \deqn{F(x) = H(u) + [1 - H(u)] G(x)}
#' where \eqn{H(x)} and \eqn{G(X)} are the BCKDE and conditional GPD
#' cumulative distribution functions respectively.
#' 
#' The cumulative distribution function for pre-specified \eqn{\phi_u}, upto the
#' threshold \eqn{x \le u}, is given by:
#' \deqn{F(x) = (1 - \phi_u) H(x)/H(u)}
#' and above the threshold \eqn{x > u}:
#' \deqn{F(x) = \phi_u + [1 - \phi_u] G(x)}
#' Notice that these definitions are equivalent when \eqn{\phi_u = 1 - H(u)}.
#' 
#' Unlike the standard KDE, there is no general rule-of-thumb bandwidth for all the
#' BCKDE, with only certain methods having a guideline in the literature, so none
#' have been implemented. Hence, a bandwidth must always be specified and you should
#' consider using \code{\link[evmix:fbckdengpd]{fbckdengpd}} of 
#' \code{\link[evmix:fbckden]{fbckden}} function for cross-validation
#' MLE for bandwidth.
#' 
#' See \code{\link[evmix:gpd]{gpd}} for details of GPD upper tail component and 
#' \code{\link[evmix:bckden]{dbckden}} for details of BCKDE bulk component.
#'
#' @section Boundary Correction Methods:
#' 
#' See \code{\link[evmix:bckden]{dbckden}} for details of BCKDE methods.
#' 
#' @section Warning:
#' The \code{"simple"}, \code{"renorm"}, \code{"beta1"}, \code{"beta2"}, \code{"gamma1"} 
#' and \code{"gamma2"} boundary correction methods may require renormalisation using
#' numerical integration which can be very slow. In particular, the numerical integration
#' is extremely slow for the \code{kernel="uniform"}, due to the adaptive quadrature in
#' the \code{\link[stats:integrate]{integrate}} function
#' being particularly slow for functions with step-like behaviour.
#' 
#' @return \code{\link[evmix:bckdengpd]{dbckdengpd}} gives the density, 
#' \code{\link[evmix:bckdengpd]{pbckdengpd}} gives the cumulative distribution function,
#' \code{\link[evmix:bckdengpd]{qbckdengpd}} gives the quantile function and 
#' \code{\link[evmix:bckdengpd]{rbckdengpd}} gives a random sample.
#' 
#' @note Unlike most of the other extreme value mixture model functions the 
#' \code{\link[evmix:bckdengpd]{bckdengpd}} functions have not been vectorised as
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
#' \code{\link[evmix:bckdengpd]{rbckdengpd}} is 1.
#' 
#' The \code{xmax} option is only relevant for the beta and copula methods, so a
#' warning is produced if this is not \code{NULL} for in other methods.
#' The \code{offset} option is only relevant for the \code{"logtrans"} method, so a
#' warning is produced if this is not \code{NULL} for in other methods.
#' 
#' Missing (\code{NA}) and Not-a-Number (\code{NaN}) values in \code{x},
#' \code{p} and \code{q} are passed through as is and infinite values are set to
#' \code{NA}. None of these are not permitted for the parameters or kernel centres.
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
#' MacDonald, A., C. J. Scarrott, and D. S. Lee (2011). Boundary correction, consistency
#' and robustness of kernel densities using extreme value theory. Submitted.
#' Available from: \url{http://www.math.canterbury.ac.nz/~c.scarrott}.
#' 
#' Wand, M. and Jones, M.C. (1995). Kernel Smoothing. Chapman && Hall.
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}. 
#' 
#' @section Acknowledgments: Based on code
#' by Anna MacDonald produced for MATLAB.
#' 
#' @seealso \code{\link[evmix:gpd]{gpd}}, \code{\link[evmix:kernels]{kernels}}, 
#' \code{\link[evmix:kfun]{kfun}},
#' \code{\link[stats:density]{density}}, \code{\link[stats:bandwidth]{bw.nrd0}}
#' and \code{\link[ks:kde.1d]{dkde}} in \code{\link[ks:kde.1d]{ks}} package.
#' 
#' @aliases bckdengpd dbckdengpd pbckdengpd qbckdengpd rbckdengpd
#' @family  kden kdengpd kdengpdcon bckden bckdengpd bckdengpdcon
#'          fkden fkdengpd fkdengpdcon fbckden fbckdengpd fbckdengpdcon
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' par(mfrow = c(2, 2))
#' 
#' kerncentres=rgamma(500, shape = 1, scale = 2)
#' xx = seq(-0.1, 10, 0.01)
#' hist(kerncentres, breaks = 100, freq = FALSE)
#' lines(xx, dbckdengpd(xx, kerncentres, lambda = 0.5, bcmethod = "reflect"),
#' xlab = "x", ylab = "f(x)")
#' abline(v = quantile(kerncentres, 0.9))
#' 
#' plot(xx, pbckdengpd(xx, kerncentres, lambda = 0.5, bcmethod = "reflect"),
#' xlab = "x", ylab = "F(x)", type = "l")
#' lines(xx, pbckdengpd(xx, kerncentres, lambda = 0.5, xi = 0.3, bcmethod = "reflect"),
#' xlab = "x", ylab = "F(x)", col = "red")
#' lines(xx, pbckdengpd(xx, kerncentres, lambda = 0.5, xi = -0.3, bcmethod = "reflect"),
#' xlab = "x", ylab = "F(x)", col = "blue")
#' legend("topleft", paste("xi =",c(0, 0.3, -0.3)),
#'       col=c("black", "red", "blue"), lty = 1, cex = 0.5)
#'
#' kerncentres = rweibull(1000, 2, 1)
#' x = rbckdengpd(1000, kerncentres, lambda = 0.1, phiu = TRUE, bcmethod = "reflect")
#' xx = seq(0.01, 3.5, 0.01)
#' hist(x, breaks = 100, freq = FALSE)         
#' lines(xx, dbckdengpd(xx, kerncentres, lambda = 0.1, phiu = TRUE, bcmethod = "reflect"),
#' xlab = "x", ylab = "f(x)")
#'
#' lines(xx, dbckdengpd(xx, kerncentres, lambda = 0.1, xi=-0.2, phiu = 0.1, bcmethod = "reflect"),
#' xlab = "x", ylab = "f(x)", col = "red")
#' lines(xx, dbckdengpd(xx, kerncentres, lambda = 0.1, xi=0.2, phiu = 0.1, bcmethod = "reflect"),
#' xlab = "x", ylab = "f(x)", col = "blue")
#' legend("topleft", c("xi = 0", "xi = 0.2", "xi = -0.2"),
#'       col=c("black", "red", "blue"), lty = 1)
#' }
#'
NULL

#' @export
#' @aliases bckdengpd dbckdengpd pbckdengpd qbckdengpd rbckdengpd
#' @rdname  bckdengpd

# probability density function for boundary corrected kernel density estimators for the bulk
# distribution upto the threshold and conditional GPD above threshold.
dbckdengpd <- function(x, kerncentres, lambda = NULL, u = as.vector(quantile(kerncentres, 0.9)),
  sigmau = sqrt(6*var(kerncentres))/pi, xi = 0, phiu = TRUE, bw  = NULL, kernel = "gaussian",
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
  check.posparam(u)
  check.posparam(sigmau)
  check.param(xi)  
  check.phiu(phiu)
  check.logic(log)

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

    if (u >= xmax) stop("threshold must be below xmax")
    
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

  check.inputn(c(length(lambda), length(u), length(sigmau), length(xi), length(phiu)), allowscalar = TRUE) # scalar only
  
  pu = pbckden(u, kerncentres, lambda, kernel = kernel, bcmethod = bcmethod,
               proper = proper, nn = nn, offset = offset, xmax = xmax)
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

  if (nb > 0) {
    d[whichb] = log(phib) + log(sapply(x[whichb], FUN = dbckden, kerncentres = kerncentres,
                                       lambda = lambda, kernel = kernel, bcmethod = bcmethod,
                                       proper = proper, nn = nn, offset = offset, xmax = xmax))
  }
  if (nu > 0) d[whichu] = log(phiu) + dgpd(x[whichu], u, sigmau, xi, log = TRUE)
  
  if (!log) d = exp(d)
  
  d
}

#' @export
#' @aliases bckdengpd dbckdengpd pbckdengpd qbckdengpd rbckdengpd
#' @rdname  bckdengpd

# cumulative distribution function for boundary corrected kernel density estimators for the bulk
# distribution upto the threshold and conditional GPD above threshold.
pbckdengpd <- function(q, kerncentres, lambda = NULL, u = as.vector(quantile(kerncentres, 0.9)),
  sigmau = sqrt(6*var(kerncentres))/pi, xi = 0, phiu = TRUE, bw  = NULL, kernel = "gaussian",
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
  check.posparam(u)
  check.posparam(sigmau)
  check.param(xi)
  check.phiu(phiu)
  check.logic(lower.tail)

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

    if (u >= xmax) stop("threshold must be below xmax")
    
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
  
  check.inputn(c(length(lambda), length(u), length(sigmau), length(xi), length(phiu)), allowscalar = TRUE) # scalar only

  pu = pbckden(u, kerncentres, lambda, kernel = kernel, bcmethod = bcmethod,
               proper = proper, nn = nn, offset = offset, xmax = xmax)
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

  if (nb > 0) {
    p[whichb] = phib*sapply(q[whichb], FUN = pbckden, kerncentres = kerncentres, lambda = lambda, 
                            kernel = kernel, bcmethod = bcmethod,
                            proper = proper, nn = nn, offset = offset, xmax = xmax)
  }
  if (nu > 0) p[whichu] = (1 - phiu) + phiu*pgpd(q[whichu], u, sigmau, xi)
  
  if (!lower.tail) p = 1 - p
  
  p
}

#' @export
#' @aliases bckdengpd dbckdengpd pbckdengpd qbckdengpd rbckdengpd
#' @rdname  bckdengpd

# inverse cumulative distribution function for boundary corrected kernel density estimators
# for the bulk distribution upto the threshold and conditional GPD above threshold.
qbckdengpd <- function(p, kerncentres, lambda = NULL, u = as.vector(quantile(kerncentres, 0.9)),
  sigmau = sqrt(6*var(kerncentres))/pi, xi = 0, phiu = TRUE, bw  = NULL, kernel = "gaussian",
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
  check.posparam(u)
  check.posparam(sigmau)
  check.param(xi)
  check.phiu(phiu)
  check.logic(lower.tail)
    
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

    if (u >= xmax) stop("threshold must be below xmax")
    
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

  check.inputn(c(length(lambda), length(u), length(sigmau), length(xi), length(phiu)), allowscalar = TRUE) # scalar only

  if (!lower.tail) p = 1 - p
  
  pu = pbckden(u, kerncentres, lambda, kernel = kernel, bcmethod = bcmethod,
               proper = proper, nn = nn, offset = offset, xmax = xmax)
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
  
  if (nb > 0) {
    q[whichb] = qbckden(p[whichb] / phib, kerncentres, lambda, kernel = kernel, 
                        bcmethod = bcmethod, proper = proper, nn = nn, offset = offset, xmax = xmax)
  }
  if (nu > 0) q[whichu] = qgpd(p[whichu], u, sigmau, xi, phiu)
  
  q
  
}

#' @export
#' @aliases bckdengpd dbckdengpd pbckdengpd qbckdengpd rbckdengpd
#' @rdname  bckdengpd

# random number generation for boundary corrected kernel density estimators
# for the bulk distribution upto the threshold and conditional GPD above threshold.
rbckdengpd <- function(n = 1, kerncentres, lambda = NULL, u = as.vector(quantile(kerncentres, 0.9)),
  sigmau = sqrt(6*var(kerncentres))/pi, xi = 0, phiu = TRUE, bw  = NULL, kernel = "gaussian",
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
  check.posparam(u)
  check.posparam(sigmau)
  check.param(xi)
  check.phiu(phiu)
  
  qbckdengpd(runif(n), kerncentres, lambda, u, sigmau, xi, phiu, bw, kernel, bcmethod,
    proper, nn, offset, xmax)
}
