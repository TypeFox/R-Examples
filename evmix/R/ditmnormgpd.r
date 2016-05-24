#' @name itmnormgpd
#' 
#' @title Normal Bulk and GPD Tail Interval Transition Mixture Model
#'
#' @description Density, cumulative distribution function, quantile function and
#'   random number generation for the normal bulk and GPD tail 
#'   interval transition mixture model. The
#'   parameters are the normal mean \code{nmean} and standard deviation \code{nsd},
#'   threshold \code{u}, interval half-width \code{epsilon}, GPD scale
#'   \code{sigmau} and shape \code{xi}.
#'
#' @param epsilon interval half-width
#' @inheritParams normgpd
#' @inheritParams gpd
#' 
#' @details The interval transition mixture model combines a normal for
#'   the bulk model with GPD for the tail model, with a smooth transition
#'   over the interval \eqn{(u-epsilon, u+epsilon)}. The mixing function warps
#'   the normal to map from \eqn{(u-epsilon, u)} to \eqn{(u-epsilon, u+epsilon)} and
#'   warps the GPD from \eqn{(u, u+epsilon)} to \eqn{(u-epsilon, u+epsilon)}.
#'   
#'   The cumulative distribution function is defined by 
#'   \deqn{F(x)=\kappa(H_t(q(x)) + G(p(x)))}
#'   where \eqn{H_t(x)} and \eqn{G(x)} are the truncated normal and
#'   conditional GPD cumulative distribution functions 
#'   (i.e. \code{pnorm(x, nmean, nsd)} and
#'    \code{pgpd(x, u, sigmau, xi)}) respectively. The truncated 
#'    normal is not renormalised to be proper, so \eqn{H_t(x)} contrubutes
#'    \code{pnorm(u, nmean, nsd)} to the cdf for all \eqn{x\geq (u + \epsilon)}.
#'    The normalisation constant \eqn{\kappa} ensures a proper density, given by 
#'    \code{1/(1+pnorm(u, nmean, nsd))} where 1 is from GPD component and
#'    latter is contribution from normal component.
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
#'   \eqn{q(x)=u} for all \eqn{x\ge u}, so the cdf/pdf of the normal model can be used
#'   directly. We do not have to define cdf/pdf for the non-proper truncated normal
#'   seperately. As such \eqn{q'(x)=0} for all \eqn{x\ge u} in
#'   \code{\link[evmix:internal]{qmixxprime}}, which also makes it clearer that
#'   normal does not contribute to the tail above the interval and vice-versa. 
#'   
#'   The quantile function within the transition interval is not available in
#'   closed form, so has to be solved numerically. Outside of the
#'   interval, the quantile are obtained from the normal and GPD components directly.
#' 
#' @return \code{\link[evmix:itmnormgpd]{ditmnormgpd}} gives the density, 
#' \code{\link[evmix:itmnormgpd]{pitmnormgpd}} gives the cumulative distribution function,
#' \code{\link[evmix:itmnormgpd]{qitmnormgpd}} gives the quantile function and 
#' \code{\link[evmix:itmnormgpd]{ritmnormgpd}} gives a random sample.
#' 
#' @note All inputs are vectorised except \code{log} and \code{lower.tail}.
#' The main inputs (\code{x}, \code{p} or \code{q}) and parameters must be either
#' a scalar or a vector. If vectors are provided they must all be of the same length,
#' and the function will be evaluated for each element of vector. In the case of 
#' \code{\link[evmix:itmnormgpd]{ritmnormgpd}} any input vector must be of length \code{n}.
#' 
#' Default values are provided for all inputs, except for the fundamentals 
#' \code{x}, \code{q} and \code{p}. The default sample size for 
#' \code{\link[evmix:itmnormgpd]{ritmnormgpd}} is 1.
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
#' Holden, L. and Haug, O. (2013). A mixture model for unsupervised tail
#' estimation. arxiv:0902.4137
#' 
#' @author Alfadino Akbar and Carl Scarrott \email{carl.scarrott@@canterbury.ac.nz}
#'
#' @seealso \code{\link[evmix:normgpd]{normgpd}}, \code{\link[evmix:gpd]{gpd}}
#'    and \code{\link[stats:Normal]{dnorm}}
#' @aliases itmnormgpd ditmnormgpd pitmnormgpd qitmnormgpd ritmnormgpd
#' @family  litmnormgpd fitmnormgpd
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' par(mfrow = c(2, 2))
#' 
#' xx = seq(-4, 5, 0.01)
#' u = 1.5
#' epsilon = 0.4
#' kappa = 1/(1 + pnorm(u, 0, 1))
#' 
#' f = ditmnormgpd(xx, nmean = 0, nsd = 1, epsilon, u, sigmau = 1, xi = 0.5)
#' plot(xx, f, ylim = c(0, 1), xlim = c(-4, 5), type = 'l', lwd = 2, xlab = "x", ylab = "density")
#' lines(xx, kappa * dgpd(xx, u, sigmau = 1, xi = 0.5), col = "red", lty = 2, lwd = 2)
#' lines(xx, kappa * dnorm(xx, 0, 1), col = "blue", lty = 2, lwd = 2)
#' abline(v = u + epsilon * seq(-1, 1), lty = c(2, 1, 2))
#' legend('topright', c('Normal-GPD ITM', 'kappa*Normal', 'kappa*GPD'),
#'       col = c("black", "blue", "red"), lty = c(1, 2, 2), lwd = 2)
#' 
#' # cdf contributions
#' F = pitmnormgpd(xx, nmean = 0, nsd = 1, epsilon, u, sigmau = 1, xi = 0.5)
#' plot(xx, F, ylim = c(0, 1), xlim = c(-4, 5), type = 'l', lwd = 2, xlab = "x", ylab = "cdf")
#' lines(xx[xx > u], kappa * (pnorm(u, 0, 1) + pgpd(xx[xx > u], u, sigmau = 1, xi = 0.5)),
#'      col = "red", lty = 2, lwd = 2)
#' lines(xx[xx <= u], kappa * pnorm(xx[xx <= u], 0, 1), col = "blue", lty = 2, lwd = 2)
#' abline(v = u + epsilon * seq(-1, 1), lty = c(2, 1, 2))
#' legend('topleft', c('Normal-GPD ITM', 'kappa*Normal', 'kappa*GPD'),
#'       col = c("black", "blue", "red"), lty = c(1, 2, 2), lwd = 2)
#'
#' # simulated data density histogram and overlay true density 
#' x = ritmnormgpd(10000, nmean = 0, nsd = 1, epsilon, u, sigmau = 1, xi = 0.5)
#' hist(x, freq = FALSE, breaks = seq(-4, 1000, 0.1), xlim = c(-4, 5))
#' lines(xx, ditmnormgpd(xx, nmean = 0, nsd = 1, epsilon, u, sigmau = 1, xi = 0.5),
#'   lwd = 2, col = 'black')  
#' }
#' 
NULL

#' @export
#' @aliases itmnormgpd ditmnormgpd pitmnormgpd qitmnormgpd ritmnormgpd
#' @rdname  itmnormgpd

# probability density function for normal bulk with GPD tail
# interval transition mixture model
ditmnormgpd = function(x, nmean = 0, nsd = 1, epsilon = nsd, u = qnorm(0.9, nmean, nsd),
                sigmau = nsd, xi = 0, log = FALSE) {
  
  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(nmean, allowvec = TRUE)
  check.posparam(nsd, allowvec = TRUE)
  check.posparam(epsilon, allowvec = TRUE, allowzero = TRUE)
  check.param(u, allowvec = TRUE)
  check.posparam(sigmau, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.logic(log)
  
  n = check.inputn(c(length(x), length(nmean), length(nsd), length(epsilon),
                     length(u), length(sigmau), length(xi)), allowscalar = TRUE)
  oneparam = (check.inputn(c(length(nmean), length(nsd), length(epsilon),
                             length(u), length(sigmau), length(xi)), allowscalar = TRUE) == 1)
  
  if (any(is.infinite(x))) warning("infinite quantiles set to NA")

  x[is.infinite(x)] = NA # user will have to deal with infinite cases
  
  x = rep(x, length.out = n)
  nmean = rep(nmean, length.out = n)  
  nsd = rep(nsd, length.out = n)
  epsilon = rep(epsilon, length.out = n)
  u = rep(u, length.out = n)
  sigmau = rep(sigmau, length.out = n)
  xi = rep(xi, length.out = n)
  
  # normalisation constant
  kappa = 1/(1 + pnorm(u, nmean, nsd))
  
  d = x # pass through NA/NaN as entered
  
  whichnonmiss = which(!is.na(x))
  
  # separate out case of scalar parameters in which this only needs to be calculated once
  if (oneparam) {
    d[whichnonmiss] = kappa[1]*(dnorm(qmix(x[whichnonmiss], u[1], epsilon[1]),
                      nmean[1], nsd[1]) * qmixprime(x[whichnonmiss], u[1], epsilon[1])
                      + dgpd(-qmix(-x[whichnonmiss], -u[1], epsilon[1]),
                      u[1], sigmau[1], xi[1]) * qmixprime(-x[whichnonmiss], -u[1], epsilon[1]))
  } else {
    nok = length(whichnonmiss)
    for (i in 1:nok) {
      d[whichnonmiss[i]] = kappa[whichnonmiss[i]]*(dnorm(qmix(x[whichnonmiss[i]], u[whichnonmiss[i]], epsilon[whichnonmiss[i]]),
                           nmean[whichnonmiss[i]], nsd[whichnonmiss[i]]) * qmixprime(x[whichnonmiss[i]], u[whichnonmiss[i]], epsilon[whichnonmiss[i]])
                           + dgpd(-qmix(-x[whichnonmiss[i]], -u[whichnonmiss[i]], epsilon[whichnonmiss[i]]),
                           u[whichnonmiss[i]], sigmau[whichnonmiss[i]], xi[whichnonmiss[i]]) * qmixprime(-x[whichnonmiss[i]], -u[whichnonmiss[i]], epsilon[whichnonmiss[i]]))
    }
  }
  if (log) d = log(d)
  
  d
}

#' @export
#' @aliases itmnormgpd ditmnormgpd pitmnormgpd qitmnormgpd ritmnormgpd
#' @rdname  itmnormgpd

# cumulative distribution function for normal bulk with GPD tail
# interval transition mixture model
pitmnormgpd = function(q, nmean = 0, nsd = 1, epsilon = nsd, u = qnorm(0.9, nmean, nsd),
                sigmau = nsd, xi = 0, lower.tail = TRUE) {
  
  # Check properties of inputs
  check.quant(q, allowna = TRUE, allowinf = TRUE)
  check.param(nmean, allowvec = TRUE)
  check.posparam(nsd, allowvec = TRUE)
  check.posparam(epsilon, allowvec = TRUE, allowzero = TRUE)
  check.param(u, allowvec = TRUE)
  check.posparam(sigmau, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.logic(lower.tail)
  
  n = check.inputn(c(length(q), length(nmean), length(nsd), length(epsilon),
                     length(u), length(sigmau), length(xi)), allowscalar = TRUE)
  oneparam = (check.inputn(c(length(nmean), length(nsd), length(epsilon),
                             length(u), length(sigmau), length(xi)), allowscalar = TRUE) == 1)
  
  if (any(is.infinite(q))) warning("infinite quantiles set to NA")
  
  q[is.infinite(q)] = NA # user will have to deal with infinite cases
  
  q = rep(q, length.out = n)
  nmean = rep(nmean, length.out = n)  
  nsd = rep(nsd, length.out = n)
  epsilon = rep(epsilon, length.out = n)
  u = rep(u, length.out = n)
  sigmau = rep(sigmau, length.out = n)
  xi = rep(xi, length.out = n)
  
  # normalisation constant
  kappa = 1/(1 + pnorm(u, nmean, nsd))

  p = q # will pass through NA/NaN as entered
  
  whichnonmiss = which(!is.na(q))
  
  # separate out case of scalar parameters in which this only needs to be calculated once
  if (oneparam) {
    p[whichnonmiss] = kappa[whichnonmiss]*(pnorm(qmix(q[whichnonmiss], u[1], epsilon[1]), nmean[1], nsd[1])
                      + pgpd(-qmix(-q[whichnonmiss], -u[1], epsilon[1]), u[1], sigmau[1], xi[1]))
  } else {
    nok = length(whichnonmiss)
    for (i in 1:nok) {
      p[whichnonmiss[i]] = kappa[whichnonmiss[i]]*(pnorm(qmix(q[whichnonmiss[i]], u[whichnonmiss[i]], epsilon[whichnonmiss[i]]), nmean[whichnonmiss[i]], nsd[whichnonmiss[i]])
                           + pgpd(-qmix(-q[whichnonmiss[i]], -u[whichnonmiss[i]], epsilon[whichnonmiss[i]]), u[whichnonmiss[i]], sigmau[whichnonmiss[i]], xi[whichnonmiss[i]]))
    }
  }
  
  if (!lower.tail) p = 1 - p
  
  p
}

#' @export
#' @aliases itmnormgpd ditmnormgpd pitmnormgpd qitmnormgpd ritmnormgpd
#' @rdname  itmnormgpd

# inverse cumulative distribution function for normal bulk with GPD tail
# interval transition mixture model
qitmnormgpd = function(p, nmean = 0, nsd = 1, epsilon = nsd, u = qnorm(0.9, nmean, nsd),
                sigmau = nsd, xi = 0, lower.tail = TRUE) {
  
  # Check properties of inputs
  check.prob(p, allowna = TRUE)
  check.param(nmean, allowvec = TRUE)
  check.posparam(nsd, allowvec = TRUE)
  check.posparam(epsilon, allowvec = TRUE, allowzero = TRUE)
  check.param(u, allowvec = TRUE)
  check.posparam(sigmau, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  check.logic(lower.tail)
  
  n = check.inputn(c(length(p), length(nmean), length(nsd), length(epsilon),
                     length(u), length(sigmau), length(xi)), allowscalar = TRUE)
  
  if (!lower.tail) p = 1 - p
    
  p = rep(p, length.out = n)
  nmean = rep(nmean, length.out = n)  
  nsd = rep(nsd, length.out = n)
  epsilon = rep(epsilon, length.out = n)
  u = rep(u, length.out = n)
  sigmau = rep(sigmau, length.out = n)
  xi = rep(xi, length.out = n)
  
  # No closed form solution for quantile function within the interval (u-epsilon, u+epsilon) exclusive,
  # need to solve numerically
  pdmmmin = function(q, cprob, nmean, nsd, epsilon, u, sigmau, xi) {
    
    cdfmm = pitmnormgpd(q, nmean, nsd, epsilon, u, sigmau, xi)
    
    if (!is.finite(cdfmm)) {
      qdiff = 1e6
    } else {
      qdiff = abs(cdfmm - cprob)
    }
    qdiff
  }
  
  findqdmm = function(cprob, interval, nmean, nsd, epsilon, u, sigmau, xi) { 
    
    gt = try(optimize(f = pdmmmin, interval = interval,
                      cprob, nmean, nsd, epsilon, u, sigmau, xi)$minimum)
    
    if (inherits(gt, "try-error")) {
      gt = NA
    }
    gt
  }
  
  # normalisation constant
  kappa = 1/(1 + pnorm(u, nmean, nsd))
  
  q = p # will pass through NA/NaN as entered
  
  whichnonmiss = which(!is.na(p))
  
  # those below and above interval (inclusive) are directly resolved
  # those in interval (exclusive) need numerical solver
  plb = kappa * pnorm(u - epsilon, nmean, nsd)
  plu =  kappa * (pnorm(u, nmean, nsd) + pgpd(u + epsilon, u, sigmau, xi))
  whichb = whichnonmiss[which(p[whichnonmiss] <= plb[whichnonmiss])]
  whichi = whichnonmiss[which((p[whichnonmiss] > plb[whichnonmiss]) & (p[whichnonmiss] < plu[whichnonmiss]))]
  whichu = whichnonmiss[which(p[whichnonmiss] >= plu[whichnonmiss])]
  nlb = length(whichb)
  nli = length(whichi)
  nlu = length(whichu)
  
  if (nlb > 0) q[whichb] = qnorm(p[whichb]/kappa[whichb], nmean[whichb], nsd[whichb])
  if (nlu > 0) q[whichu] = qgpd(p[whichu]/kappa[whichu] - pnorm(u[whichu], nmean[whichu], nsd[whichu]), u[whichu], sigmau[whichu], xi[whichu])

  if (nli > 0) {
    for (i in 1:length(whichi)) {
      q[whichi[i]] = findqdmm(p[whichi[i]], u[whichi[i]] + epsilon[whichi[i]] * c(-1, 1),
                    nmean[whichi[i]], nsd[whichi[i]], epsilon[whichi[i]], 
                    u[whichi[i]], sigmau[whichi[i]], xi[whichi[i]])
    }
  }
  
  q                     
}

#' @export
#' @aliases itmnormgpd ditmnormgpd pitmnormgpd qitmnormgpd ritmnormgpd
#' @rdname  itmnormgpd

# random number generation for normal bulk with GPD tail
# interval transition mixture model
ritmnormgpd = function(n = 1, nmean = 0, nsd = 1, epsilon = nsd, u = qnorm(0.9, nmean, nsd),
                sigmau = nsd, xi = 0) {
  
  # Check properties of inputs
  check.n(n)
  check.param(nmean, allowvec = TRUE)
  check.posparam(nsd, allowvec = TRUE)
  check.posparam(epsilon, allowvec = TRUE, allowzero = TRUE)
  check.param(u, allowvec = TRUE)
  check.posparam(sigmau, allowvec = TRUE)
  check.param(xi, allowvec = TRUE)
  
  check.inputn(c(n, length(nmean), length(nsd), length(epsilon), 
                 length(u), length(sigmau), length(xi)), allowscalar = TRUE)
  
  if (any(xi == 1)) stop("shape cannot be 1")
  
  qitmnormgpd(runif(n), nmean, nsd, epsilon, u, sigmau, xi)
}

