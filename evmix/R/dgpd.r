#' @name gpd
#' 
#' @title Generalised Pareto Distribution (GPD)
#'
#' @description Density, cumulative distribution function, quantile function and
#'   random number generation for the generalised Pareto distribution, either
#'   as a conditional on being above the threshold \code{u} or unconditional. 
#'
#' @param x          quantiles
#' @param q          quantiles
#' @param p          cumulative probabilities
#' @param n          sample size (positive integer)
#' @param u          threshold
#' @param sigmau     scale parameter (positive)
#' @param xi         shape parameter
#' @param phiu       probability of being above threshold \eqn{[0, 1]}
#' @param log        logical, if TRUE then log density
#' @param lower.tail logical, if FALSE then upper tail probabilities
#' 
#' @details The GPD with parameters scale \eqn{\sigma_u} and shape \eqn{\xi} has
#' conditional density of being above the threshold \code{u} given by
#'
#'  \deqn{f(x | X > u) = 1/\sigma_u [1 + \xi(x - u)/\sigma_u]^{-1/\xi - 1}}
#'  
#' for non-zero \eqn{\xi}, \eqn{x > u} and \eqn{\sigma_u > 0}. Further, 
#' \eqn{[1+\xi (x - u) / \sigma_u] > 0} which for \eqn{\xi < 0} implies 
#' \eqn{u < x \le u - \sigma_u/\xi}. In the special case of \eqn{\xi = 0}
#' considered in the limit \eqn{\xi \rightarrow 0}, which is
#' treated here as \eqn{|\xi| < 1e-6}, it reduces to the exponential:
#' 
#' \deqn{f(x | X > u) = 1/\sigma_u exp(-(x - u)/\sigma_u).}
#' 
#' The unconditional density is obtained by mutltiplying this by the
#' survival probability (or \emph{tail fraction}) \eqn{\phi_u = P(X > u)}
#' giving \eqn{f(x) = \phi_u f(x | X > u)}.
#' 
#' The syntax of these functions are similar to those of the 
#' \code{\link[evd:gpd]{evd}} package, so most code using these functions can
#' be reused. The key difference is the introduction of \code{phiu} to
#' permit output of unconditional quantities.
#' 
#' @return \code{\link[evmix:gpd]{dgpd}} gives the density,
#' \code{\link[evmix:gpd]{pgpd}} gives the cumulative distribution function,
#' \code{\link[evmix:gpd]{qgpd}} gives the quantile function and 
#' \code{\link[evmix:gpd]{rgpd}} gives a random sample.
#' 
#' @note All inputs are vectorised except \code{log} and \code{lower.tail}.
#' The main inputs (\code{x}, \code{p} or \code{q}) and parameters must be either
#' a scalar or a vector. If vectors are provided they must all be of the same length,
#' and the function will be evaluated for each element of vector. In the case of 
#' \code{\link[evmix:gpd]{rgpd}} any input vector must be of length \code{n}.
#' 
#' Default values are provided for all inputs, except for the fundamentals 
#' \code{x}, \code{q} and \code{p}. The default threshold \code{u=0} and tail fraction
#' \code{phiu=1} which essentially assumes the user provide excesses above 
#' \code{u} by default, rather than exceedances. The default sample size for 
#' \code{\link[evmix:gpd]{rgpd}} is 1.
#' 
#' Missing (\code{NA}) and Not-a-Number (\code{NaN}) values in \code{x},
#' \code{p} and \code{q} are passed through as is and infinite values are set to
#' \code{NA}. None of these are not permitted for the parameters.
#' 
#' Some key differences arise for \code{phiu=1} and \code{phiu<1} (see examples below):
#' 
#' \enumerate{
#' \item For \code{phiu=1} the \code{\link[evmix:gpd]{dgpd}} evaluates as zero for
#' quantiles below the threshold \code{u} and \code{\link[evmix:gpd]{pgpd}}
#' evaluates over \eqn{[0, 1]}.
#' 
#' \item For \code{phiu=1} then \code{\link[evmix:gpd]{pgpd}} evaluates as zero
#' below the threshold \code{u}. For \code{phiu<1} it evaluates as \eqn{1-\phi_u} at
#' the threshold and \code{NA} below the threshold.
#' 
#' \item For \code{phiu=1} the quantiles from \code{\link[evmix:gpd]{qgpd}} are
#' above threshold and equal to threshold for \code{phiu=0}. For \code{phiu<1} then
#' within upper tail, \code{p > 1 - phiu}, it will give conditional quantiles
#' above threshold, but when below the threshold, \code{p <= 1 - phiu}, these
#' are set to \code{NA}.
#' 
#' \item When simulating GPD variates using \code{\link[evmix:gpd]{rgpd}} if
#' \code{phiu=1} then all values are above the threshold. For \code{phiu<1} then
#' a standard uniform \eqn{U} is simulated and the variate will be classified as
#' above the threshold if \eqn{u<\phi}, and below the threshold otherwise. This is
#' equivalent to a binomial random variable for simulated number of exceedances. Those
#' above the threshold are then simulated from the conditional GPD and those below
#' the threshold and set to \code{NA}.
#' }
#' 
#' These conditions are intuitive and consistent with \code{\link[evd:gpd]{evd}},
#' which assumes missing data are below threshold.
#' 
#' Error checking of the inputs (e.g. invalid probabilities) is carried out and
#' will either stop or give warning message as appropriate.
#' 
#' @references
#' 
#' \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
#' 
#' Coles, S.G. (2001). An Introduction to Statistical Modelling of Extreme Values.
#' Springer Series in Statistics. Springer-Verlag: London.
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @section Acknowledgments: Based on the
#' \code{\link[evd:gpd]{gpd}} functions in the \code{\link[evd:gpd]{evd}} package for which their author's contributions are gratefully acknowledged.
#' They are designed to have similar syntax and functionality to simplify the transition for users of these packages.
#'   
#' @seealso \code{\link[evd:gpd]{evd}} package and \code{\link[evd:fpot]{fpot}}
#' 
#' @aliases gpd dgpd pgpd qgpd rgpd
#' @family  gpd fgpd
#' 
#' @examples
#' set.seed(1)
#' par(mfrow = c(2, 2))
#' 
#' x = rgpd(1000) # simulate sample from GPD
#' xx = seq(-1, 10, 0.01)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-1, 10))
#' lines(xx, dgpd(xx))
#'
#' # three tail behaviours
#' plot(xx, pgpd(xx), type = "l")
#' lines(xx, pgpd(xx, xi = 0.3), col = "red")
#' lines(xx, pgpd(xx, xi = -0.3), col = "blue")
#' legend("bottomright", paste("xi =",c(0, 0.3, -0.3)),
#'   col=c("black", "red", "blue"), lty = 1)
#' 
#' # GPD when xi=0 is exponential, and demonstrating phiu
#' x = rexp(1000)
#' hist(x, breaks = 100, freq = FALSE, xlim = c(-1, 10))
#' lines(xx, dgpd(xx, u = 0, sigmau = 1, xi = 0), lwd = 2)
#' lines(xx, dgpd(xx, u = 0.5, phiu = 1 - pexp(0.5)), col = "red", lwd = 2)
#' lines(xx, dgpd(xx, u = 1.5, phiu = 1 - pexp(1.5)), col = "blue", lwd = 2)
#' legend("topright", paste("u =",c(0, 0.5, 1.5)),
#'   col=c("black", "red", "blue"), lty = 1, lwd = 2)
#' 
#' # Quantile function and phiu
#' p = pgpd(xx)
#' plot(qgpd(p), p, type = "l")
#' lines(xx, pgpd(xx, u = 2), col = "red")
#' lines(xx, pgpd(xx, u = 5, phiu = 0.2), col = "blue")
#' legend("bottomright", c("u = 0 phiu = 1","u = 2 phiu = 1","u = 5 phiu = 0.2"),
#'   col=c("black", "red", "blue"), lty = 1)
#'   
NULL

#' @export
#' @aliases gpd dgpd pgpd qgpd rgpd
#' @rdname  gpd

# probability density function for GPD
dgpd <- function(x, u = 0, sigmau = 1, xi = 0, phiu = 1, log = FALSE) {

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(u, allowvec = TRUE)
  check.posparam(sigmau, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.prob(phiu) # don't use check.phiu as TRUE only valid for mixture models
  check.logic(log)

  n = check.inputn(c(length(x), length(u), length(sigmau), length(xi), length(phiu)), allowscalar = TRUE)

  if (any(is.infinite(x))) warning("infinite quantiles set to NA")

  x[is.infinite(x)] = NA # user will have to deal with infinite cases
 
  x = rep(x, length.out = n)
  u = rep(u, length.out = n)
  sigmau = rep(sigmau, length.out = n)
  xi = rep(xi, length.out = n)
  phiu = rep(phiu, length.out = n)
  
  d = x # will pass through NA/NaN as entered
  
  yu = (x - u)/sigmau # used when shape is zero
  syu = 1 + xi*yu     # used when shape non-zero
  
  # check for x values in range
  yind = ((yu > 0) & (syu > 0))

  d[which(!yind)] = log(0) # zero density is default
  
  # special case when xi parameter is zero (or close to it)
  shape0ind = abs(xi) < 1e-6
  nshape0 = sum(shape0ind)

  if (nshape0 > 0) {
    whichexp = which(shape0ind & yind)
    d[whichexp] = -log(sigmau[whichexp]) - yu[whichexp]
  }
  if (nshape0 < n) {
    whichxi = which(!shape0ind & yind)
    d[whichxi] = -log(sigmau[whichxi]) - (1/xi[whichxi] + 1) * log(syu[whichxi])
  }
  
  d = d + log(phiu) # unconditional density
  
  if (!log) d = exp(d)
  
  d
}

#' @export
#' @aliases gpd dgpd pgpd qgpd rgpd
#' @rdname  gpd

# cumulative distribution function for GPD
pgpd <- function(q, u = 0, sigmau = 1, xi = 0, phiu = 1, lower.tail = TRUE) {

  # Check properties of inputs
  check.quant(q, allowna = TRUE, allowinf = TRUE)
  check.param(u, allowvec = TRUE)
  check.posparam(sigmau, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.prob(phiu) # don't use check.phiu as TRUE only valid for mixture models
  check.logic(lower.tail)

  n = check.inputn(c(length(q), length(u), length(sigmau), length(xi), length(phiu)), allowscalar = TRUE)

  if (any(is.infinite(q))) warning("infinite quantiles set to NA")

  q[is.infinite(q)] = NA # user will have to deal with infinite cases
 
  q = rep(q, length.out = n)
  u = rep(u, length.out = n)
  sigmau = rep(sigmau, length.out = n)
  xi = rep(xi, length.out = n)
  phiu = rep(phiu, length.out = n)
  
  yu = pmax(q - u, 0) / sigmau # used when shape is zero
  syu = pmax(1 + xi*yu, 0)     # used when shape non-zero
  
  # check for x values in range
  yind = (yu > 0)

  p = q # will pass through NA/NaN as entered
  
  p[which(!yind)] = ifelse(phiu[which(!yind)] == 1, 0, NA) # NA is default, but if only GPD then gives 0 below threshold
  
  # special case when xi parameter is zero (or close to it)
  shape0ind = abs(xi) < 1e-6
  nshape0 = sum(shape0ind)
    
  if (nshape0 > 0) {
    whichexp = which(shape0ind & yind)
    p[whichexp] = 1 - exp(-yu[whichexp])
  }
  if (nshape0 < n) {
    whichxi = which(!shape0ind & yind)
    p[whichxi] = 1 - syu[whichxi] ^ (-1 / xi[whichxi])
  }
  
  p = 1 - (1 - p) * phiu
  
  if (!lower.tail) p = 1 - p
  
  p
}

#' @export
#' @aliases gpd dgpd pgpd qgpd rgpd
#' @rdname  gpd

# inverse cumulative distribution function for GPD
qgpd <- function(p, u = 0, sigmau = 1, xi = 0, phiu = 1, lower.tail = TRUE) {

  # Check properties of inputs
  check.prob(p, allowna = TRUE)
  check.param(u, allowvec = TRUE)
  check.posparam(sigmau, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.prob(phiu) # don't use check.phiu as TRUE only valid for mixture models
  check.logic(lower.tail)

  n = check.inputn(c(length(p), length(u), length(sigmau), length(xi), length(phiu)), allowscalar = TRUE)
 
  if (!lower.tail) p = 1 - p

  p = rep(p, length.out = n)
  u = rep(u, length.out = n)
  sigmau = rep(sigmau, length.out = n)
  xi = rep(xi, length.out = n)
  phiu = rep(phiu, length.out = n)
  
  # check for x values in range
  yind = (p > (1 - phiu)) 

  q = p # will pass through NA/NaN as entered
  
  q[which(!yind)] = ifelse(phiu[which(!yind)] == 1, u[which(!yind)], NA) # NA is default, but if only GPD then gives threshold

  # special case when xi parameter is zero (or close to it)
  shape0ind = abs(xi) < 1e-6
  nshape0 = sum(shape0ind)
    
  if (nshape0 > 0) {
    whichexp = which(shape0ind & yind)
    q[whichexp] = u[whichexp] - sigmau[whichexp] * log((1 - p[whichexp]) / phiu[whichexp])
  }
  if (nshape0 < n) {
    whichxi = which(!shape0ind & yind)
    q[whichxi] = u[whichxi] + sigmau[whichxi] * (((1 - p[whichxi]) / phiu[whichxi]) ^ (-xi[whichxi]) - 1) / xi[whichxi]
  }
    
  q
}

#' @export
#' @aliases gpd dgpd pgpd qgpd rgpd
#' @rdname  gpd

# random number generation for GPD
rgpd <- function(n = 1, u = 0, sigmau = 1, xi = 0, phiu = 1) {

  # Check properties of inputs
  check.n(n)
  check.param(u, allowvec = TRUE)
  check.posparam(sigmau, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.prob(phiu) # don't use check.phiu as TRUE only valid for mixture models

  n = check.inputn(c(n, length(u), length(sigmau), length(xi), length(phiu)), allowscalar = TRUE)

  if (any(xi == 1)) stop("shape cannot be 1")
    
  u = rep(u, length.out = n)
  sigmau = rep(sigmau, length.out = n)
  xi = rep(xi, length.out = n)
  phiu = rep(phiu, length.out = n)

  yind = (runif(n) < phiu)
  
  # special case when xi parameter is zero (or close to it)
  shape0ind = abs(xi) < 1e-6
  nshape0 = sum(shape0ind)
  
  r = rep(NA, n)
  
  if (nshape0 > 0) {
    whichexp = which(shape0ind & yind)
    r[whichexp] = u[whichexp] + sigmau[whichexp] * rexp(u[whichexp])
  }
  if (nshape0 < n) {
    whichxi = which(!shape0ind  & yind)
    r[whichxi] = u[whichxi] + sigmau[whichxi] * (runif(length(whichxi)) ^ (-xi[whichxi]) - 1) / xi[whichxi]
  }
    
  r
}
