#' @name itmweibullgpd
#' 
#' @title Weibull Bulk and GPD Tail Interval Transition Mixture Model
#'
#' @description Density, cumulative distribution function, quantile function and
#'   random number generation for the Weibull bulk and GPD tail 
#'   interval transition mixture model. The
#'   parameters are the Weibull shape \code{wshape} and scale \code{wscale},
#'   threshold \code{u}, interval half-width \code{epsilon}, GPD scale
#'   \code{sigmau} and shape \code{xi}.
#'
#' @inheritParams itmnormgpd
#' @inheritParams weibullgpd
#' @inheritParams gpd
#' 
#' @details The interval transition mixture model combines a Weibull for
#'   the bulk model with GPD for the tail model, with a smooth transition
#'   over the interval \eqn{(u-epsilon, u+epsilon)}. The mixing function warps
#'   the Weibull to map from \eqn{(u-epsilon, u)} to \eqn{(u-epsilon, u+epsilon)} and
#'   warps the GPD from \eqn{(u, u+epsilon)} to \eqn{(u-epsilon, u+epsilon)}.
#'   
#'   The cumulative distribution function is defined by 
#'   \deqn{F(x)=\kappa(H_t(q(x)) + G(p(x)))}
#'   where \eqn{H_t(x)} and \eqn{G(X)} are the truncated Weibull and
#'   conditional GPD cumulative distribution functions 
#'   (i.e. \code{pweibull(x, wshape, wscale)} and
#'    \code{pgpd(x, u, sigmau, xi)}) respectively. The truncated 
#'    Weibull is not renormalised to be proper, so \eqn{H_t(x)} contrubutes
#'    \code{pweibull(u, wshape, wscale)} to the cdf for all \eqn{x\geq (u + \epsilon)}.
#'    The normalisation constant \eqn{\kappa} ensures a proper density, given by 
#'    \code{1/(1+pweibull(u, wshape, wscale))} where 1 is from GPD component and
#'    latter is contribution from Weibull component.
#'   
#'   The mixing functions \eqn{q(x)} and \eqn{p(x)} suggested by Holden and Haug (2013)
#'   have been implemented. These are symmetric about the threshold \eqn{u}. So for
#'   computational convenience only \eqn{q(x;u)} has been implemented as 
#'   \code{\link[evmix:internal]{qmix}}
#'   for a given \eqn{u}, with the complementary mixing function is then defined as
#'   \eqn{p(x;u)=-q(-x;-u)}.
#'   
#'   A minor adaptation of the mixing function has been applied.  For the mixture model to
#'   function correctly \eqn{q(x)>=u} for all \eqn{x\ge u+\epsilon}, as then the bulk model will contribute
#'   the constant \eqn{H_t(u)=H(u)} for all \eqn{x} above the interval. Holden and Haug (2013) define
#'   \eqn{q(x)=x-\epsilon} for all \eqn{x\ge u}. For more straightforward and interpretable 
#'   computational implementation the mixing function has been set to the threshold
#'   \eqn{q(x)=u} for all \eqn{x\ge u}, so the cdf/pdf of the Weibull model can be used
#'   directly. We do not have to define cdf/pdf for the non-proper truncated Weibull
#'   seperately. As such \eqn{q'(x)=0} for all \eqn{x\ge u} in
#'   \code{\link[evmix:internal]{qmixxprime}}, which also it makes clearer that
#'   Weibull does not contribute to the tail above the interval and vice-versa. 
#'   
#'   The quantile function within the transition interval is not available in
#'   closed form, so has to be solved numerically. Outside of the
#'   interval, the quantile are obtained from the Weibull and GPD components directly.
#' 
#' @return \code{\link[evmix:itmweibullgpd]{ditmweibullgpd}} gives the density, 
#' \code{\link[evmix:itmweibullgpd]{pitmweibullgpd}} gives the cumulative distribution function,
#' \code{\link[evmix:itmweibullgpd]{qitmweibullgpd}} gives the quantile function and 
#' \code{\link[evmix:itmweibullgpd]{ritmweibullgpd}} gives a random sample.
#' 
#' @note All inputs are vectorised except \code{log} and \code{lower.tail}.
#' The main inputs (\code{x}, \code{p} or \code{q}) and parameters must be either
#' a scalar or a vector. If vectors are provided they must all be of the same length,
#' and the function will be evaluated for each element of vector. In the case of 
#' \code{\link[evmix:itmweibullgpd]{ritmweibullgpd}} any input vector must be of length \code{n}.
#' 
#' Default values are provided for all inputs, except for the fundamentals 
#' \code{x}, \code{q} and \code{p}. The default sample size for 
#' \code{\link[evmix:itmweibullgpd]{ritmweibullgpd}} is 1.
#' 
#' Missing (\code{NA}) and Not-a-Number (\code{NaN}) values in \code{x},
#' \code{p} and \code{q} are passed through as is and infinite values are set to
#' \code{NA}. None of these are not permitted for the parameters.
#' 
#' Error checking of the inputs (e.g. invalid probabilities) is carried out and
#' will either stop or give warning message as appropriate.
#' 
#' @references
#' \url{http://en.wikipedia.org/wiki/Weibull_distribution}
#' 
#' \url{http://en.wikipedia.org/wiki/Generalized_Pareto_distribution}
#' 
#' Scarrott, C.J. and MacDonald, A. (2012). A review of extreme value
#' threshold estimation and uncertainty quantification. REVSTAT - Statistical
#' Journal 10(1), 33-59. Available from \url{http://www.ine.pt/revstat/pdf/rs120102.pdf}
#' 
#' Holden, L. and Haug, O. (2013). A mixture model for unsupervised tail
#' estimation. arxiv:0902.4137
#' 
#' @author Alfadino Akbar and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @seealso \code{\link[evmix:weibullgpd]{weibullgpd}}, \code{\link[evmix:gpd]{gpd}}
#'    and \code{\link[stats:Weibull]{dweibull}}
#' @aliases itmweibullgpd ditmweibullgpd pitmweibullgpd qitmweibullgpd ritmweibullgpd
#' @family  litmweibullgpd fitmweibullgpd
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' par(mfrow = c(2, 2))
#' 
#' xx = seq(0.001, 5, 0.01)
#' u = 1.5
#' epsilon = 0.4
#' kappa = 1/(1 + pweibull(u, 2, 1))
#' 
#' f = ditmweibullgpd(xx, wshape = 2, wscale = 1, epsilon, u, sigmau = 1, xi = 0.5)
#' plot(xx, f, ylim = c(0, 1), xlim = c(0, 5), type = 'l', lwd = 2, xlab = "x", ylab = "density")
#' lines(xx, kappa * dgpd(xx, u, sigmau = 1, xi = 0.5), col = "red", lty = 2, lwd = 2)
#' lines(xx, kappa * dweibull(xx, 2, 1), col = "blue", lty = 2, lwd = 2)
#' abline(v = u + epsilon * seq(-1, 1), lty = c(2, 1, 2))
#' legend('topright', c('Weibull-GPD ITM', 'kappa*Weibull', 'kappa*GPD'),
#'       col = c("black", "blue", "red"), lty = c(1, 2, 2), lwd = 2)
#' 
#' # cdf contributions
#' F = pitmweibullgpd(xx, wshape = 2, wscale = 1, epsilon, u, sigmau = 1, xi = 0.5)
#' plot(xx, F, ylim = c(0, 1), xlim = c(0, 5), type = 'l', lwd = 2, xlab = "x", ylab = "cdf")
#' lines(xx[xx > u], kappa * (pweibull(u, 2, 1) + pgpd(xx[xx > u], u, sigmau = 1, xi = 0.5)),
#'      col = "red", lty = 2, lwd = 2)
#' lines(xx[xx <= u], kappa * pweibull(xx[xx <= u], 2, 1), col = "blue", lty = 2, lwd = 2)
#' abline(v = u + epsilon * seq(-1, 1), lty = c(2, 1, 2))
#' legend('topright', c('Weibull-GPD ITM', 'kappa*Weibull', 'kappa*GPD'),
#'       col = c("black", "blue", "red"), lty = c(1, 2, 2), lwd = 2)
#'
#' # simulated data density histogram and overlay true density 
#' x = ritmweibullgpd(10000, wshape = 2, wscale = 1, epsilon, u, sigmau = 1, xi = 0.5)
#' hist(x, freq = FALSE, breaks = seq(0, 1000, 0.1), xlim = c(0, 5))
#' lines(xx, ditmweibullgpd(xx, wshape = 2, wscale = 1, epsilon, u, sigmau = 1, xi = 0.5),
#'   lwd = 2, col = 'black')  
#' }
#' 
NULL

#' @export
#' @aliases itmweibullgpd ditmweibullgpd pitmweibullgpd qitmweibullgpd ritmweibullgpd
#' @rdname  itmweibullgpd

# probability density function for Weibull bulk with GPD tail
# interval transition mixture model
ditmweibullgpd = function(x, wshape = 1, wscale = 1, 
                epsilon = sqrt(wscale^2 * gamma(1 + 2/wshape) - (wscale * gamma(1 + 1/wshape))^2),          
                u = qweibull(0.9, wshape, wscale),
                sigmau = sqrt(wscale^2 * gamma(1 + 2/wshape) - (wscale * gamma(1 + 1/wshape))^2),
                xi = 0, log = FALSE) {
  
  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.posparam(wshape, allowvec = TRUE)
  check.posparam(wscale, allowvec = TRUE)
  check.posparam(epsilon, allowvec = TRUE, allowzero = TRUE)
  check.posparam(u, allowvec = TRUE)
  check.posparam(sigmau, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.logic(log)
  
  n = check.inputn(c(length(x), length(wshape), length(wscale), length(epsilon),
                     length(u), length(sigmau), length(xi)), allowscalar = TRUE)
  oneparam = (check.inputn(c(length(wshape), length(wscale), length(epsilon),
                             length(u), length(sigmau), length(xi)), allowscalar = TRUE) == 1)
  
  if (any(is.infinite(x))) warning("infinite quantiles set to NA")

  x[is.infinite(x)] = NA # user will have to deal with infinite cases
  
  x = rep(x, length.out = n)
  wshape = rep(wshape, length.out = n)  
  wscale = rep(wscale, length.out = n)
  epsilon = rep(epsilon, length.out = n)
  u = rep(u, length.out = n)
  sigmau = rep(sigmau, length.out = n)
  xi = rep(xi, length.out = n)
  
  # normalisation constant
  kappa = 1/(1 + pweibull(u, wshape, wscale))
  
  d = x # pass through NA/NaN as entered
  
  whichnonmiss = which(!is.na(x))
  
  # separate out case of scalar parameters in which this only needs to be calculated once
  if (oneparam) {
    d[whichnonmiss] = kappa[1]*(dweibull(qmix(x[whichnonmiss], u[1], epsilon[1]),
                      wshape[1], wscale[1]) * qmixprime(x[whichnonmiss], u[1], epsilon[1])
                      + dgpd(-qmix(-x[whichnonmiss], -u[1], epsilon[1]),
                      u[1], sigmau[1], xi[1]) * qmixprime(-x[whichnonmiss], -u[1], epsilon[1]))
  } else {
    nok = length(whichnonmiss)
    for (i in 1:nok) {
      d[whichnonmiss[i]] = kappa[whichnonmiss[i]]*(dweibull(qmix(x[whichnonmiss[i]], u[whichnonmiss[i]], epsilon[whichnonmiss[i]]),
                           wshape[whichnonmiss[i]], wscale[whichnonmiss[i]]) * qmixprime(x[whichnonmiss[i]], u[whichnonmiss[i]], epsilon[whichnonmiss[i]])
                           + dgpd(-qmix(-x[whichnonmiss[i]], -u[whichnonmiss[i]], epsilon[whichnonmiss[i]]),
                           u[whichnonmiss[i]], sigmau[whichnonmiss[i]], xi[whichnonmiss[i]]) * qmixprime(-x[whichnonmiss[i]], -u[whichnonmiss[i]], epsilon[whichnonmiss[i]]))
    }
  }
  if (log) d = log(d)
  
  d
}

#' @export
#' @aliases itmweibullgpd ditmweibullgpd pitmweibullgpd qitmweibullgpd ritmweibullgpd
#' @rdname  itmweibullgpd

# cumulative distribution function for Weibull bulk with GPD tail
# interval transition mixture model
pitmweibullgpd = function(q, wshape = 1, wscale = 1, 
                epsilon = sqrt(wscale^2 * gamma(1 + 2/wshape) - (wscale * gamma(1 + 1/wshape))^2),
                u = qweibull(0.9, wshape, wscale),
                sigmau = sqrt(wscale^2 * gamma(1 + 2/wshape) - (wscale * gamma(1 + 1/wshape))^2),
                xi = 0, lower.tail = TRUE) {
  
  # Check properties of inputs
  check.quant(q, allowna = TRUE, allowinf = TRUE)
  check.posparam(wshape, allowvec = TRUE)
  check.posparam(wscale, allowvec = TRUE)
  check.posparam(epsilon, allowvec = TRUE, allowzero = TRUE)
  check.posparam(u, allowvec = TRUE)
  check.posparam(sigmau, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.logic(lower.tail)
  
  n = check.inputn(c(length(q), length(wshape), length(wscale), length(epsilon),
                     length(u), length(sigmau), length(xi)), allowscalar = TRUE)
  oneparam = (check.inputn(c(length(wshape), length(wscale), length(epsilon),
                             length(u), length(sigmau), length(xi)), allowscalar = TRUE) == 1)
  
  if (any(is.infinite(q))) warning("infinite quantiles set to NA")
  
  q[is.infinite(q)] = NA # user will have to deal with infinite cases
  
  q = rep(q, length.out = n)
  wshape = rep(wshape, length.out = n)  
  wscale = rep(wscale, length.out = n)
  epsilon = rep(epsilon, length.out = n)
  u = rep(u, length.out = n)
  sigmau = rep(sigmau, length.out = n)
  xi = rep(xi, length.out = n)
  
  # normalisation constant
  kappa = 1/(1 + pweibull(u, wshape, wscale))

  p = q # will pass through NA/NaN as entered
  
  whichnonmiss = which(!is.na(q))
  
  # separate out case of scalar parameters in which this only needs to be calculated once
  if (oneparam) {
    p[whichnonmiss] = kappa[whichnonmiss]*(pweibull(qmix(q[whichnonmiss], u[1], epsilon[1]), wshape[1], wscale[1])
                      + pgpd(-qmix(-q[whichnonmiss], -u[1], epsilon[1]), u[1], sigmau[1], xi[1]))
  } else {
    nok = length(whichnonmiss)
    for (i in 1:nok) {
      p[whichnonmiss[i]] = kappa[whichnonmiss[i]]*(pweibull(qmix(q[whichnonmiss[i]], u[whichnonmiss[i]], epsilon[whichnonmiss[i]]), wshape[whichnonmiss[i]], wscale[whichnonmiss[i]])
                           + pgpd(-qmix(-q[whichnonmiss[i]], -u[whichnonmiss[i]], epsilon[whichnonmiss[i]]), u[whichnonmiss[i]], sigmau[whichnonmiss[i]], xi[whichnonmiss[i]]))
    }
  }
  
  if (!lower.tail) p = 1 - p
  
  p
}

#' @export
#' @aliases itmweibullgpd ditmweibullgpd pitmweibullgpd qitmweibullgpd ritmweibullgpd
#' @rdname  itmweibullgpd

# inverse cumulative distribution function for Weibull bulk with GPD tail
# interval transition mixture model
qitmweibullgpd = function(p, wshape = 1, wscale = 1, 
                epsilon = sqrt(wscale^2 * gamma(1 + 2/wshape) - (wscale * gamma(1 + 1/wshape))^2),
                u = qweibull(0.9, wshape, wscale),
                sigmau = sqrt(wscale^2 * gamma(1 + 2/wshape) - (wscale * gamma(1 + 1/wshape))^2),
                xi = 0, lower.tail = TRUE) {
  
  # Check properties of inputs
  check.prob(p, allowna = TRUE)
  check.posparam(wshape, allowvec = TRUE)
  check.posparam(wscale, allowvec = TRUE)
  check.posparam(epsilon, allowvec = TRUE, allowzero = TRUE)
  check.posparam(u, allowvec = TRUE)
  check.posparam(sigmau, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.logic(lower.tail)
  
  n = check.inputn(c(length(p), length(wshape), length(wscale), length(epsilon),
                     length(u), length(sigmau), length(xi)), allowscalar = TRUE)
  
  if (!lower.tail) p = 1 - p
    
  p = rep(p, length.out = n)
  wshape = rep(wshape, length.out = n)  
  wscale = rep(wscale, length.out = n)
  epsilon = rep(epsilon, length.out = n)
  u = rep(u, length.out = n)
  sigmau = rep(sigmau, length.out = n)
  xi = rep(xi, length.out = n)
  
  # No closed form solution for quantile function within the interval (u-epsilon, u+epsilon) exclusive,
  # need to solve numerically
  pdmmmin = function(q, cprob, wshape, wscale, epsilon, u, sigmau, xi) {
    
    cdfmm = pitmweibullgpd(q, wshape, wscale, epsilon, u, sigmau, xi)
    
    if (!is.finite(cdfmm)) {
      qdiff = 1e6
    } else {
      qdiff = abs(cdfmm - cprob)
    }
    qdiff
  }
  
  findqdmm = function(cprob, interval, wshape, wscale, epsilon, u, sigmau, xi) {
    
    gt = try(optimize(f = pdmmmin, interval = interval,
                      cprob, wshape, wscale, epsilon, u, sigmau, xi)$minimum)
    
    if (inherits(gt, "try-error")) {
      gt = NA
    }
    gt
  }
  
  # normalisation constant
  kappa = 1/(1 + pweibull(u, wshape, wscale))
  
  q = p # will pass through NA/NaN as entered
  
  whichnonmiss = which(!is.na(p))
  
  # those below and above interval (inclusive) are directly resolved
  # those in interval (exclusive) need numerical solver
  plb = kappa * pweibull(u - epsilon, wshape, wscale)
  plu =  kappa * (pweibull(u, wshape, wscale) + pgpd(u + epsilon, u, sigmau, xi))
  whichb = whichnonmiss[which(p[whichnonmiss] <= plb[whichnonmiss])]
  whichi = whichnonmiss[which((p[whichnonmiss] > plb[whichnonmiss]) & (p[whichnonmiss] < plu[whichnonmiss]))]
  whichu = whichnonmiss[which(p[whichnonmiss] >= plu[whichnonmiss])]
  nlb = length(whichb)
  nli = length(whichi)
  nlu = length(whichu)
  
  if (nlb > 0) q[whichb] = qweibull(p[whichb]/kappa[whichb], wshape[whichb], wscale[whichb])
  if (nlu > 0) q[whichu] = qgpd(p[whichu]/kappa[whichu] - pweibull(u[whichu], wshape[whichu], wscale[whichu]), u[whichu], sigmau[whichu], xi[whichu])

  if (nli > 0) {
    for (i in 1:length(whichi)) {    
      q[whichi[i]] = findqdmm(p[whichi[i]], u[whichi[i]] + epsilon[whichi[i]] * c(-1, 1),
                    wshape[whichi[i]], wscale[whichi[i]], epsilon[whichi[i]], 
                    u[whichi[i]], sigmau[whichi[i]], xi[whichi[i]])
    }
  }
  
  q                     
}

#' @export
#' @aliases itmweibullgpd ditmweibullgpd pitmweibullgpd qitmweibullgpd ritmweibullgpd
#' @rdname  itmweibullgpd

# random number generation for Weibull bulk with GPD tail
# interval transition mixture model
ritmweibullgpd = function(n = 1, wshape = 1, wscale = 1, 
                epsilon = sqrt(wscale^2 * gamma(1 + 2/wshape) - (wscale * gamma(1 + 1/wshape))^2),
                u = qweibull(0.9, wshape, wscale),
                sigmau = sqrt(wscale^2 * gamma(1 + 2/wshape) - (wscale * gamma(1 + 1/wshape))^2),
                xi = 0) {
  
  # Check properties of inputs
  check.n(n)
  check.posparam(wshape, allowvec = TRUE)
  check.posparam(wscale, allowvec = TRUE)
  check.posparam(epsilon, allowvec = TRUE, allowzero = TRUE)
  check.posparam(u, allowvec = TRUE)
  check.posparam(sigmau, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  
  check.inputn(c(n, length(wshape), length(wscale), length(epsilon), 
                 length(u), length(sigmau), length(xi)), allowscalar = TRUE)
  
  if (any(xi == 1)) stop("shape cannot be 1")
  
  qitmweibullgpd(runif(n), wshape, wscale, epsilon, u, sigmau, xi)
}

