#' @name hpd
#' 
#' @title Hybrid Pareto Extreme Value Mixture Model
#'
#' @description Density, cumulative distribution function, quantile function and
#'   random number generation for the hybrid Pareto extreme value mixture model.
#'   The parameters are the normal mean \code{nmean} and standard deviation \code{nsd} and 
#'   GPD shape \code{xi}.  
#'
#' @inheritParams normgpd
#' @inheritParams gpd
#' 
#' @details Extreme value mixture model combining normal distribution for the bulk
#' below the threshold and GPD for upper tail which is continuous in its zeroth and
#' first derivative at the threshold. 
#' 
#' But it has one important difference to all the other mixture models. The
#' hybrid Pareto does not include the usual tail fraction \code{phiu} scaling, 
#' i.e. so the GPD is not treated as a conditional model for the exceedances. 
#' The unscaled GPD is simply spliced with the normal truncated at the
#' threshold, with no rescaling to account for the proportion above the
#' threshold being applied. The parameters have to adjust for the lack of tail 
#' fraction scaling.
#' 
#' The cumulative distribution function defined upto the 
#' threshold \eqn{x \le u}, given by:
#' \deqn{F(x) = H(x) / r }
#' and above the threshold \eqn{x > u}:
#' \deqn{F(x) = (H(u) +  G(x)) / r }
#' where \eqn{H(x)} and \eqn{G(X)} are the normal and conditional GPD
#' cumulative distribution functions. The normalisation constant \eqn{r} ensures a proper
#' density and is given by\code{r = 1 + pnorm(u, mean = nmean, sd = nsd)}, i.e. the 1 comes from
#' integration of the unscaled GPD and the second term is from the usual normal component.
#' 
#' The two continuity constraints leads to the threshold \code{u} and GPD scale \code{sigmau} being replaced
#' by a function of the normal mean, standard deviation and GPD shape parameters. 
#' Determined from setting \eqn{h(u) = g(u)} where \eqn{h(x)} and \eqn{g(x)} are the normal and unscaled GPD
#' density functions (i.e. \code{dnorm(u, nmean, nsd)} and
#' \code{dgpd(u, u, sigmau, xi)}). The continuity constraint on its first derivative at the threshold 
#' means that \eqn{h'(u) = g'(u)}. Then the Lambert-W function is used for replacing
#' the threshold u and GPD scale sigmau in terms of the normal mean, standard deviation
#' and GPD shape xi.
#' 
#' See \code{\link[evmix:gpd]{gpd}} for details of GPD upper tail component and 
#'\code{\link[stats:Normal]{dnorm}} for details of normal bulk component.
#' 
#' @return \code{\link[evmix:hpd]{dhpd}} gives the density, 
#' \code{\link[evmix:hpd]{phpd}} gives the cumulative distribution function,
#' \code{\link[evmix:hpd]{qhpd}} gives the quantile function and 
#' \code{\link[evmix:hpd]{rhpd}} gives a random sample.
#' 
#' @note All inputs are vectorised except \code{log} and \code{lower.tail}.
#' The main inputs (\code{x}, \code{p} or \code{q}) and parameters must be either
#' a scalar or a vector. If vectors are provided they must all be of the same length,
#' and the function will be evaluated for each element of vector. In the case of 
#' \code{\link[evmix:hpd]{rhpd}} any input vector must be of length \code{n}.
#' 
#' Default values are provided for all inputs, except for the fundamentals 
#' \code{x}, \code{q} and \code{p}. The default sample size for 
#' \code{\link[evmix:hpd]{rhpd}} is 1.
#' 
#' Missing (\code{NA}) and Not-a-Number (\code{NaN}) values in \code{x},
#' \code{p} and \code{q} are passed through as is and infinite values are set to
#' \code{NA}. None of these are not permitted for the parameters.
#' 
#' Error checking of the inputs (e.g. invalid probabilities) is carried out and
#' will either stop or give warning message as appropriate.
#' 
#' @references
#' \url{http://en.wikipedia.org/wiki/Normal_distribution}
#' 
#' \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
#' 
#' Scarrott, C.J. and MacDonald, A. (2012). A review of extreme value
#' threshold estimation and uncertainty quantification. REVSTAT - Statistical
#' Journal 10(1), 33-59. Available from \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}
#' 
#' Carreau, J. and Y. Bengio (2008). A hybrid Pareto model for asymmetric fat-tailed data:
#' the univariate case. Extremes 12 (1), 53-76.
#' 
#' @author Yang Hu and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @seealso \code{\link[evmix:gpd]{gpd}} and \code{\link[stats:Normal]{dnorm}}.
#' 
#' The \code{\link[condmixt:condmixt-package]{condmixt}} package written by one of the
#' original authors of the hybrid Pareto model (Carreau and Bengio, 2008) also has 
#' similar functions for the hybrid Pareto \code{\link[condmixt:hpareto]{hpareto}} and
#' mixture of hybrid Paretos \code{\link[condmixt:hparetomixt]{hparetomixt}}, which are
#' more flexible as they also permit the model to be truncated at zero.
#' 
#' @aliases hpd dhpd phpd qhpd rhpd
#' @family  hpd hpdcon fhpd fhpdcon normgpd normgpdcon fnormgpd fnormgpdcon
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' par(mfrow = c(2, 2))
#' 
#' xx = seq(-5, 20, 0.01)
#' f1 = dhpd(xx, nmean = 0, nsd = 1, xi = 0.4)
#' plot(xx, f1, type = "l")
#' abline(v = 0.4942921)
#' 
#' # three tail behaviours
#' plot(xx, phpd(xx), type = "l")
#' lines(xx, phpd(xx, xi = 0.3), col = "red")
#' lines(xx, phpd(xx, xi = -0.3), col = "blue")
#' legend("bottomright", paste("xi =",c(0, 0.3, -0.3)),
#'   col=c("black", "red", "blue"), lty = 1)
#'  
#' sim = rhpd(10000, nmean = 0, nsd = 1.5, xi = 0.2)
#' hist(sim, freq = FALSE, 100, xlim = c(-5, 20), ylim = c(0, 0.2))
#' lines(xx, dhpd(xx, nmean = 0, nsd = 1.5, xi = 0.2), col = "blue")
#' 
#' plot(xx, dhpd(xx, nmean = 0, nsd = 1.5, xi = 0), type = "l")
#' lines(xx, dhpd(xx, nmean = 0, nsd = 1.5, xi = 0.2), col = "red")
#' lines(xx, dhpd(xx, nmean = 0, nsd = 1.5, xi = -0.2), col = "blue")
#' legend("topright", c("xi = 0", "xi = 0.2", "xi = -0.2"),
#'   col=c("black", "red", "blue"), lty = 1)
#' }
#' 
NULL

#' @export
#' @aliases hpd dhpd phpd qhpd rhpd
#' @rdname  hpd

# probability density function for hybrid Pareto model
dhpd <- function(x, nmean = 0, nsd = 1, xi = 0, log = FALSE) {
  
  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(nmean, allowvec = TRUE)
  check.posparam(nsd, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.logic(log)

  n = check.inputn(c(length(x), length(nmean), length(nsd), length(xi)), allowscalar = TRUE)

  if (any(is.infinite(x))) warning("infinite quantiles set to NA")

  x[is.infinite(x)] = NA # user will have to deal with infinite cases
      
  x = rep(x, length.out = n)
  nmean = rep(nmean, length.out = n)
  nsd = rep(nsd, length.out = n)
  xi = rep(xi, length.out = n)
  
  z = (1 + xi)^2/(2*pi)
  wz = lambert_W0(z)
  
  u = nmean + nsd * sqrt(wz) * sign(1 + xi)
  sigmau = nsd * abs(1 + xi) / sqrt(wz)
  
  check.posparam(sigmau, allowvec = TRUE)
  check.param(u, allowvec = TRUE)
  
  r = 1 + pnorm(u, nmean, nsd)
  
  d = x # will pass through NA/NaN as entered

  whichb = which(x <= u)
  nb = length(whichb)
  whichu = which(x > u)
  nu = length(whichu)
  
  if (nb > 0) d[whichb] = dnorm(x[whichb], nmean[whichb], nsd[whichb], log = TRUE) - log(r[whichb])
  if (nu > 0) d[whichu] = dgpd(x[whichu], u[whichu], sigmau[whichu], xi[whichu], log = TRUE) - log(r[whichu])
  
  if (!log) d = exp(d)
  
  d
}

#' @export
#' @aliases hpd dhpd phpd qhpd rhpd
#' @rdname  hpd

# cumulative distribution function for hybrid Pareto model
phpd <- function(q, nmean = 0, nsd = 1, xi = 0, lower.tail = TRUE) {
  
  # Check properties of inputs
  check.quant(q, allowna = TRUE, allowinf = TRUE)
  check.param(nmean, allowvec = TRUE)
  check.posparam(nsd, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.logic(lower.tail)

  n = check.inputn(c(length(q), length(nmean), length(nsd), length(xi)), allowscalar = TRUE)

  if (any(is.infinite(q))) warning("infinite quantiles set to NA")

  q[is.infinite(q)] = NA # user will have to deal with infinite cases
    
  q = rep(q, length.out = n)
  nmean = rep(nmean, length.out = n)
  nsd = rep(nsd, length.out = n)
  xi = rep(xi, length.out = n)
  
  z = (1 + xi)^2/(2*pi)
  wz = lambert_W0(z)
  
  u = nmean + nsd * sqrt(wz) * sign(1 + xi)
  sigmau = nsd * abs(1 + xi) / sqrt(wz)
  
  check.posparam(sigmau, allowvec = TRUE)
  check.param(u, allowvec = TRUE)

  r = 1 + pnorm(u, nmean, nsd)
    
  p = q # will pass through NA/NaN as entered
  
  whichb = which(q <= u)
  nb = length(whichb)
  whichu = which(q > u)
  nu = length(whichu)
  
  if (nb > 0) p[whichb] = pnorm(q[whichb], nmean[whichb], nsd[whichb]) / r[whichb]
  if (nu > 0) p[whichu] = (pnorm(u[whichu], nmean[whichu], nsd[whichu]) + pgpd(q[whichu], u[whichu], sigmau[whichu], xi[whichu])) / r[whichu]
  
  if (!lower.tail) p = 1 - p
  
  p
}

#' @export
#' @aliases hpd dhpd phpd qhpd rhpd
#' @rdname  hpd

# inverse cumulative distribution function for hybrid Pareto model
qhpd <- function(p, nmean = 0, nsd = 1, xi = 0, lower.tail = TRUE) {
  
  # Check properties of inputs
  check.prob(p, allowna = TRUE)
  check.param(nmean, allowvec = TRUE)
  check.posparam(nsd, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.logic(lower.tail)

  n = check.inputn(c(length(p), length(nmean), length(nsd), length(xi)), allowscalar = TRUE)

  if (!lower.tail) p = 1 - p
  
  p = rep(p, length.out = n)
  nmean = rep(nmean, length.out = n)
  nsd = rep(nsd, length.out = n)
  xi = rep(xi, length.out = n)
  
  z = (1 + xi)^2/(2*pi)
  wz = lambert_W0(z)
  
  u = nmean + nsd * sqrt(wz) * sign(1 + xi)
  sigmau = nsd * abs(1 + xi) / sqrt(wz)
  
  check.posparam(sigmau, allowvec = TRUE)
  check.param(u, allowvec = TRUE)
    
  r = 1 + pnorm(u, nmean, nsd)
  
  q = p # will pass through NA/NaN as entered
  
  phi = pnorm(u, nmean, nsd)
  phiu = phi/(1 + phi)
  
  whichb = which(q <= phiu)
  nb = length(whichb)
  whichu = which(q > phiu)
  nu = length(whichu)
  
  if (nb>0) {q[whichb] = qnorm((1 + phi[whichb]) * p[whichb], nmean[whichb], nsd[whichb])}
  if (nu>0) {q[whichu] = qgpd((1 + phi[whichu]) * p[whichu] - phi[whichu], 0, sigmau[whichu], xi[whichu]) + u[whichu]}
    
  q  
}

#' @export
#' @aliases hpd dhpd phpd qhpd rhpd
#' @rdname  hpd

# random number generation for hybrid Pareto model
rhpd <- function(n = 1, nmean = 0, nsd = 1, xi = 0) {
  
  # Check properties of inputs
  check.n(n)
  check.param(nmean, allowvec = TRUE)
  check.posparam(nsd, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)

  n = check.inputn(c(n, length(nmean), length(nsd), length(xi)), allowscalar = TRUE)

  if (any(xi == 1)) stop("shape cannot be 1")
    
  qhpd(runif(n), nmean, nsd, xi)
}
