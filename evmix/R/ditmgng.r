#' @name itmgng
#' 
#' @title Normal Bulk with GPD Upper and Lower Tails Interval Transition Mixture Model

#' @description Density, cumulative distribution function, quantile function and
#'   random number generation for the extreme value mixture model with normal
#'   for bulk distribution between the upper and lower thresholds with
#'   conditional GPD's for the two tails and interval transition. The parameters are the normal mean
#'   \code{nmean} and standard deviation \code{nsd}, interval half-width \code{espilon},
#'   lower tail (threshold \code{ul}, GPD scale \code{sigmaul} and shape \code{xil} and
#'   tail fraction \code{phiul}) and upper tail (threshold \code{ur}, GPD scale
#'   \code{sigmaur} and shape \code{xiR} and tail fraction \code{phiuR}).
#'
#' @inheritParams itmnormgpd
#' @inheritParams gng
#' @inheritParams gpd
#' 
#' @details The interval transition extreme value mixture model combines a normal
#' distribution for the bulk between the lower and upper thresholds and GPD for
#' upper and lower tails, with a smooth transition over the interval 
#' \eqn{(u-epsilon, u+epsilon)} (where \eqn{u} can be exchanged for the lower and
#' upper thresholds). The mixing function warps the normal to map from 
#' \eqn{(u-epsilon, u)} to \eqn{(u-epsilon, u+epsilon)} and warps the GPD from 
#' \eqn{(u, u+epsilon)} to \eqn{(u-epsilon, u+epsilon)}.
#' 
#'   The cumulative distribution function is defined by 
#'   \deqn{F(x)=\kappa(G_l(q(x)) + H_t(r(x)) + G_u(p(x)))}
#'   where \eqn{H_t(x)} is the truncated normal cdf, i.e. \code{pnorm(x, nmean, nsd)}.
#'   The conditional GPD for the upper tail has cdf \eqn{G_u(x)}, 
#'   i.e. \code{pgpd(x, ur, sigmaur, xir)} and lower tail cdf \eqn{G_l(x)} is for the 
#'   negated support, i.e. \code{1 - pgpd(-x, -ul, sigmaul, xil)}. The truncated 
#'   normal is not renormalised to be proper, so \eqn{H_t(x)} contributes
#'   \code{pnorm(ur, nmean, nsd) - pnorm(ul, nmean, nsd)} to the cdf
#'   for all \eqn{x\geq (u_r + \epsilon)} and zero below \eqn{x\leq (u_l - \epsilon)}.
#'   The normalisation constant \eqn{\kappa} ensures a proper density, given by 
#'    \code{1/(2 + pnorm(ur, nmean, nsd) - pnorm(ul, nmean, nsd)} where the
#'    2 is from two GPD components and latter is contribution from normal component.
#'   
#'   The mixing functions \eqn{q(x)}, \eqn{r(x)} and \eqn{p(x)} are reformulated from the 
#'   \eqn{q_i(x)} suggested by Holden and Haug (2013). These are symmetric about each
#'   threshold, which for convenience will be referred to a simply \eqn{u}. So for
#'   computational convenience only a single \eqn{q(x;u)} has been implemented for the
#'   lower and upper GPD components called
#'   \code{\link[evmix:internal]{qmix}} for a given \eqn{u}, with the complementary
#'   mixing function then defined as \eqn{p(x;u)=-q(-x;-u)}. The bulk model mixing
#'   function \eqn{r(x)} utilises the equivalent of the \eqn{q(x)} for the lower threshold and
#'   \eqn{p(x)} for the upper threshold, so these are reused in the bulk mixing function  
#'   \code{\link[evmix:internal]{qgbgmix}}.
#'   
#'   A minor adaptation of the mixing function has been applied following a similar
#'   approach to that explained in \code{\link[evmix:itmnormgpd]{ditmnormgpd}}. For the
#'   bulk model mixing function \eqn{r(x)}, we need \eqn{r(x)<=ul} for all \eqn{x\le ul - epsilon} and 
#'   \eqn{r(x)>=ur} for all \eqn{x\ge ur+epsilon}, as then the bulk model will contribute
#'   zero below the lower interval and the constant \eqn{H_t(ur)=H(ur)-H(ul)} for all
#'   \eqn{x} above the upper interval. Holden and Haug (2013) define
#'   \eqn{r(x)=x-\epsilon} for all \eqn{x\ge ur} and \eqn{r(x)=x+\epsilon} for all \eqn{x\le ul}.
#'   For more straightforward and interpretable 
#'   computational implementation the mixing function has been set to the lower threshold
#'   \eqn{r(x)=u_l} for all \eqn{x\le u_l-\epsilon} and to the upper threshold
#'   \eqn{r(x)=u_r} for all \eqn{x\le u_r+\epsilon}, so the cdf/pdf of the normal model can be used
#'   directly. We do not have to define cdf/pdf for the non-proper truncated normal
#'   seperately. As such \eqn{r'(x)=0} for all \eqn{x\le u_l-\epsilon} and \eqn{x\ge u_r+\epsilon} in
#'   \code{\link[evmix:internal]{qmixxprime}}, which also makes it clearer that
#'   normal does not contribute to either tails beyond the intervals and vice-versa. 
#'   
#'   The quantile function within the transition interval is not available in
#'   closed form, so has to be solved numerically. Outside of the
#'   interval, the quantile are obtained from the normal and GPD components directly.
#'   
#' @return \code{\link[evmix:itmgng]{ditmgng}} gives the density, 
#' \code{\link[evmix:itmgng]{pitmgng}} gives the cumulative distribution function,
#' \code{\link[evmix:itmgng]{qitmgng}} gives the quantile function and 
#' \code{\link[evmix:itmgng]{ritmgng}} gives a random sample.
#' 
#' @note All inputs are vectorised except \code{log} and \code{lower.tail}.
#' The main input (\code{x}, \code{p} or \code{q}) and parameters must be either
#' a scalar or a vector. If vectors are provided they must all be of the same length,
#' and the function will be evaluated for each element of vector. In the case of 
#' \code{\link[evmix:itmgng]{ritmgng}} any input vector must be of length \code{n}.
#' 
#' Default values are provided for all inputs, except for the fundamentals 
#' \code{x}, \code{q} and \code{p}. The default sample size for 
#' \code{\link[evmix:itmgng]{ritmgng}} is 1.
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
#' @seealso   \code{\link[evmix:gng]{gng}}, \code{\link[evmix:normgpd]{normgpd}},
#'            \code{\link[evmix:gpd]{gpd}} and \code{\link[stats:Normal]{dnorm}}
#' @aliases itmgng ditmgng pitmgng qitmgng ritmgng
#' @family  normgpd normgpdcon gng gngcon fnormgpd fnormgpdcon fgng fgngcon
#' 
#' @examples
#' \dontrun{
#' set.seed(1)
#' par(mfrow = c(2, 2))
#' 
#' xx = seq(-5, 5, 0.01)
#' ul = -1.5;ur = 2
#' epsilon = 0.8
#' kappa = 1/(2 + pnorm(ur, 0, 1) - pnorm(ul, 0, 1))
#' 
#' f = ditmgng(xx, nmean = 0, nsd = 1, epsilon, ul, sigmaul = 1, xil = 0.5, ur, sigmaur = 1, xir = 0.5)
#' plot(xx, f, ylim = c(0, 0.5), xlim = c(-5, 5), type = 'l', lwd = 2, xlab = "x", ylab = "density")
#' lines(xx, kappa * dgpd(-xx, -ul, sigmau = 1, xi = 0.5), col = "blue", lty = 2, lwd = 2)
#' lines(xx, kappa * dnorm(xx, 0, 1), col = "red", lty = 2, lwd = 2)
#' lines(xx, kappa * dgpd(xx, ur, sigmau = 1, xi = 0.5), col = "green", lty = 2, lwd = 2)
#' abline(v = ul + epsilon * seq(-1, 1), lty = c(2, 1, 2), col = "blue")
#' abline(v = ur + epsilon * seq(-1, 1), lty = c(2, 1, 2), col = "green")
#' legend('topright', c('Normal-GPD ITM', 'kappa*GPD Lower', 'kappa*Normal', 'kappa*GPD Upper'),
#'       col = c("black", "blue", "red", "green"), lty = c(1, 2, 2, 2), lwd = 2)
#' 
#' # cdf contributions
#' F = pitmgng(xx, nmean = 0, nsd = 1, epsilon, ul, sigmaul = 1, xil = 0.5, ur, sigmaur = 1, xir = 0.5)
#' plot(xx, F, ylim = c(0, 1), xlim = c(-5, 5), type = 'l', lwd = 2, xlab = "x", ylab = "cdf")
#' lines(xx[xx < ul], kappa * (1 - pgpd(-xx[xx < ul], -ul, 1, 0.5)), col = "blue", lty = 2, lwd = 2)
#' lines(xx[(xx >= ul) & (xx <= ur)], kappa * (1 + pnorm(xx[(xx >= ul) & (xx <= ur)], 0, 1) -
#'       pnorm(ul, 0, 1)), col = "red", lty = 2, lwd = 2)
#' lines(xx[xx > ur], kappa * (1 + (pnorm(ur, 0, 1) - pnorm(ul, 0, 1)) +
#'       pgpd(xx[xx > ur], ur, sigmau = 1, xi = 0.5)), col = "green", lty = 2, lwd = 2)
#' abline(v = ul + epsilon * seq(-1, 1), lty = c(2, 1, 2), col = "blue")
#' abline(v = ur + epsilon * seq(-1, 1), lty = c(2, 1, 2), col = "green")
#' legend('topleft', c('Normal-GPD ITM', 'kappa*GPD Lower', 'kappa*Normal', 'kappa*GPD Upper'),
#'       col = c("black", "blue", "red", "green"), lty = c(1, 2, 2, 2), lwd = 2)
#'
#' # simulated data density histogram and overlay true density 
#' x = ritmgng(10000, nmean = 0, nsd = 1, epsilon, ul, sigmaul = 1, xil = 0.5,
#'                                                 ur, sigmaur = 1, xir = 0.5)
#' hist(x, freq = FALSE, breaks = seq(-1000, 1000, 0.1), xlim = c(-5, 5))
#' lines(xx, ditmgng(xx, nmean = 0, nsd = 1, epsilon, ul, sigmaul = 1, xil = 0.5,
#'   ur, sigmaur = 1, xir = 0.5), lwd = 2, col = 'black')
#' }
#' 
NULL

#' @export
#' @aliases itmgng ditmgng pitmgng qitmgng ritmgng
#' @rdname  itmgng

# probability density function for normal bulk with GPD's for upper and lower tails
# interval transition mixture model
ditmgng <- function(x, nmean = 0, nsd = 1, epsilon = nsd,
  ul = qnorm(0.1, nmean, nsd), sigmaul = nsd, xil = 0,
  ur = qnorm(0.9, nmean, nsd), sigmaur = nsd, xir = 0, log = FALSE) {

  # Check properties of inputs
  check.quant(x, allowna = TRUE, allowinf = TRUE)
  check.param(nmean, allowvec = TRUE)
  check.posparam(nsd, allowvec = TRUE)
  check.posparam(epsilon, allowvec = TRUE, allowzero = TRUE)
  check.param(ul, allowvec = TRUE)
  check.posparam(sigmaul, allowvec = TRUE)
  check.param(xil, allowvec = TRUE)
  check.param(ur, allowvec = TRUE)
  check.posparam(sigmaur, allowvec = TRUE)
  check.param(xir, allowvec = TRUE)
  check.logic(log)

  n = check.inputn(c(length(x), length(nmean), length(nsd), length(epsilon),
    length(ul), length(sigmaul), length(xil), length(ur), length(sigmaur), length(xir)), allowscalar = TRUE)
  oneparam = (check.inputn(c(length(nmean), length(nsd), length(epsilon),
    length(ul), length(sigmaul), length(xil), length(ur), length(sigmaur), length(xir)), allowscalar = TRUE) == 1)

  if (any(is.infinite(x))) warning("infinite quantiles set to NA")

  x[is.infinite(x)] = NA # user will have to deal with infinite cases

  if (any((ul + epsilon) >= (ur - epsilon))) stop("intervals must not overlap")
  
  x = rep(x, length.out = n)
  nmean = rep(nmean, length.out = n)
  nsd = rep(nsd, length.out = n)
  epsilon = rep(epsilon, length.out = n)
  ul = rep(ul, length.out = n)
  sigmaul = rep(sigmaul, length.out = n)
  xil = rep(xil, length.out = n)
  ur = rep(ur, length.out = n)
  sigmaur = rep(sigmaur, length.out = n)
  xir = rep(xir, length.out = n)
  
  # normalisation constant
  kappa = 1/(2 + pnorm(ur, nmean, nsd) - pnorm(ul, nmean, nsd))

  d = x # will pass through NA/NaN as entered
  
  whichnonmiss = which(!is.na(x))
  
  # separate out case of scalar parameters in which this only needs to be calculated once
  if (oneparam) {
    d[whichnonmiss] = kappa[1]*(
      dgpd(-qmix(x[whichnonmiss], ul[1], epsilon[1]), -ul[1], sigmaul[1], xil[1]) * qmixprime(x[whichnonmiss], ul[1], epsilon[1]) +
      dnorm(qgbgmix(x[whichnonmiss], ul[1], ur[1], epsilon[1]), nmean[1], nsd[1]) * qgbgmixprime(x[whichnonmiss], ul[1], ur[1], epsilon[1]) +
      dgpd(-qmix(-x[whichnonmiss], -ur[1], epsilon[1]), ur[1], sigmaur[1], xir[1]) * qmixprime(-x[whichnonmiss], -ur[1], epsilon[1]))
  } else {
    nok = length(whichnonmiss)
    for (i in 1:nok) {
      d[whichnonmiss[i]] = kappa[whichnonmiss[i]]*(
        dgpd(-qmix(x[whichnonmiss[i]], ul[whichnonmiss[i]], epsilon[whichnonmiss[i]]), -ul[whichnonmiss[i]], sigmaul[whichnonmiss[i]], xil[whichnonmiss[i]]) *
          qmixprime(x[whichnonmiss[i]], ul[whichnonmiss[i]], epsilon[whichnonmiss[i]]) +
        dnorm(qgbgmix(x[whichnonmiss[i]], ul[whichnonmiss[i]], ur[whichnonmiss[i]], epsilon[whichnonmiss[i]]), nmean[whichnonmiss[i]], nsd[whichnonmiss[i]]) *
          qgbgmixprime(x[whichnonmiss[i]], ul[whichnonmiss[i]], ur[whichnonmiss[i]], epsilon[whichnonmiss[i]]) +
        dgpd(-qmix(-x[whichnonmiss[i]], -ur[whichnonmiss[i]], epsilon[whichnonmiss[i]]), ur[whichnonmiss[i]], sigmaur[whichnonmiss[i]], xir[whichnonmiss[i]]) *
          qmixprime(-x[whichnonmiss[i]], -ur[whichnonmiss[i]], epsilon[whichnonmiss[i]]))

    }
  }
  if (log) d = log(d)

  d
}

#' @export
#' @aliases itmgng ditmgng pitmgng qitmgng ritmgng
#' @rdname  itmgng

# cumulative distribution function for normal bulk with GPD's for upper and lower tails
# interval transition mixture model
pitmgng <- function(q, nmean = 0, nsd = 1, epsilon = nsd,
  ul = qnorm(0.1, nmean, nsd), sigmaul = nsd, xil = 0,
  ur = qnorm(0.9, nmean, nsd), sigmaur = nsd, xir = 0, lower.tail = TRUE) {

  # Check properties of inputs
  check.quant(q, allowna = TRUE, allowinf = TRUE)
  check.param(nmean, allowvec = TRUE)
  check.posparam(nsd, allowvec = TRUE)
  check.posparam(epsilon, allowvec = TRUE, allowzero = TRUE)
  check.param(ul, allowvec = TRUE)
  check.posparam(sigmaul, allowvec = TRUE)
  check.param(xil, allowvec = TRUE)
  check.param(ur, allowvec = TRUE)
  check.posparam(sigmaur, allowvec = TRUE)
  check.param(xir, allowvec = TRUE)
  check.logic(lower.tail)

  n = check.inputn(c(length(q), length(nmean), length(nsd), length(epsilon),
    length(ul), length(sigmaul), length(xil), length(ur), length(sigmaur), length(xir)), allowscalar = TRUE)
  oneparam = (check.inputn(c(length(nmean), length(nsd), length(epsilon),
    length(ul), length(sigmaul), length(xil), length(ur), length(sigmaur), length(xir)), allowscalar = TRUE) == 1)

  if (any(is.infinite(q))) warning("infinite quantiles set to NA")

  q[is.infinite(q)] = NA # user will have to deal with infinite cases

  if (any((ul + epsilon) >= (ur - epsilon))) stop("intervals must not overlap")

  q = rep(q, length.out = n)
  nmean = rep(nmean, length.out = n)
  nsd = rep(nsd, length.out = n)
  epsilon = rep(epsilon, length.out = n)
  ul = rep(ul, length.out = n)
  sigmaul = rep(sigmaul, length.out = n)
  xil = rep(xil, length.out = n)
  ur = rep(ur, length.out = n)
  sigmaur = rep(sigmaur, length.out = n)
  xir = rep(xir, length.out = n)
  
  # normalisation constant
  kappa = 1/(2 + pnorm(ur, nmean, nsd) - pnorm(ul, nmean, nsd))
    
  p = q # will pass through NA/NaN as entered
  
  whichnonmiss = which(!is.na(q))
  
  # separate out case of scalar parameters in which this only needs to be calculated once
  if (oneparam) {
    p[whichnonmiss] = kappa[whichnonmiss]*(
                       (1 - pgpd(-qmix(q[whichnonmiss], ul[1], epsilon[1]), -ul[1], sigmaul[1], xil[1])) +
                       (pnorm(qgbgmix(q[whichnonmiss], ul[1], ur[1], epsilon[1]), nmean[1], nsd[1]) - 
                          pnorm(ul[1], nmean[1], nsd[1])) +
                       pgpd(-qmix(-q[whichnonmiss], -ur[1], epsilon[1]), ur[1], sigmaur[1], xir[1]))
  } else {
    nok = length(whichnonmiss)
    for (i in 1:nok) {
      p[whichnonmiss[i]] = kappa[whichnonmiss[i]]*(
          (1 - pgpd(-qmix(q[whichnonmiss[i]], ul[i], epsilon[i]), -ul[i], sigmaul[i], xil[i])) +
          (pnorm(qgbgmix(q[whichnonmiss[i]], ul[whichnonmiss[i]], ur[whichnonmiss[i]], epsilon[whichnonmiss[i]]), nmean[whichnonmiss[i]], nsd[whichnonmiss[i]]) -
             pnorm(ul[whichnonmiss[i]], nmean[whichnonmiss[i]], nsd[whichnonmiss[i]])) +
          pgpd(-qmix(-q[whichnonmiss[i]], -ur[whichnonmiss[i]], epsilon[whichnonmiss[i]]), ur[whichnonmiss[i]], sigmaur[whichnonmiss[i]], xir[whichnonmiss[i]]))
    }
  }
  
  if (!lower.tail) p = 1 - p

  p
}

#' @export
#' @aliases itmgng ditmgng pitmgng qitmgng ritmgng
#' @rdname  itmgng

# inverse cumulative distribution function for normal bulk with GPD's for upper and lower tails
# interval transition mixture model
qitmgng <- function(p, nmean = 0, nsd = 1, epsilon,
  ul = qnorm(0.1, nmean, nsd), sigmaul = nsd, xil = 0, 
  ur = qnorm(0.9, nmean, nsd), sigmaur = nsd, xir = 0, lower.tail = TRUE) {

  # Check properties of inputs
  check.prob(p, allowna = TRUE)
  check.param(nmean, allowvec = TRUE)
  check.posparam(nsd, allowvec = TRUE)
  check.posparam(epsilon, allowvec = TRUE, allowzero = TRUE)
  check.param(ul, allowvec = TRUE)
  check.posparam(sigmaul, allowvec = TRUE)
  check.param(xil, allowvec = TRUE)
  check.param(ur, allowvec = TRUE)
  check.posparam(sigmaur, allowvec = TRUE)
  check.param(xir, allowvec = TRUE)
  check.logic(lower.tail)

  n = check.inputn(c(length(p), length(nmean), length(nsd), length(epsilon),
    length(ul), length(sigmaul), length(xil), length(ur), length(sigmaur), length(xir)), allowscalar = TRUE)
  oneparam = (check.inputn(c(length(nmean), length(nsd), length(epsilon),
    length(ul), length(sigmaul), length(xil), length(ur), length(sigmaur), length(xir)), allowscalar = TRUE) == 1)

  if (any((ul + epsilon) >= (ur - epsilon))) stop("intervals must not overlap")

  if (!lower.tail) p = 1 - p

  p = rep(p, length.out = n)
  nmean = rep(nmean, length.out = n)
  nsd = rep(nsd, length.out = n)
  epsilon = rep(epsilon, length.out = n)
  ul = rep(ul, length.out = n)
  sigmaul = rep(sigmaul, length.out = n)
  xil = rep(xil, length.out = n)
  ur = rep(ur, length.out = n)
  sigmaur = rep(sigmaur, length.out = n)
  xir = rep(xir, length.out = n)
  
  
  # No closed form solution for quantile function within either of the intervals
  # (ul-epsilon, ul+epsilon) and (ur-epsilon, ur+epsilon) exclusive,
  # need to solve numerically
  pdmmmin = function(q, cprob, nmean, nsd, epsilon, ul, sigmaul, xil, ur, sigmaur, xir) {
    
    cdfmm = pitmgng(q, nmean, nsd, epsilon, ul, sigmaul, xil, ur, sigmaur, xir)
    
    if (!is.finite(cdfmm)) {
      qdiff = 1e6
    } else {
      qdiff = abs(cdfmm - cprob)
    }
    qdiff
  }
  
  findqdmm = function(cprob, interval, nmean, nsd, epsilon, ul, sigmaul, xil, ur, sigmaur, xir) {

    gt = try(optimize(f = pdmmmin, interval = interval,
                      cprob, nmean, nsd, epsilon, ul, sigmaul, xil, ur, sigmaur, xir)$minimum)
    
    if (inherits(gt, "try-error")) {
      gt = NA
    }
    gt
  }
  
  # normalisation constant
  kappa = 1/(2 + pnorm(ur, nmean, nsd) - pnorm(ul, nmean, nsd))
  
  q = p # will pass through NA/NaN as entered
  
  whichnonmiss = which(!is.na(p))
  
  # outside of both intervals (inclusive) are directly resolved
  # those in both intervals (exclusive) need numerical solver
  plul =  kappa * (1 - pgpd(-(ul - epsilon), -ul, sigmaul, xil))
  plbl = kappa * (1 + pnorm(ul + epsilon, nmean, nsd) - pnorm(ul, nmean, nsd))
  plbr = kappa * (1 + pnorm(ur - epsilon, nmean, nsd) - pnorm(ul, nmean, nsd))
  plur =  kappa * (1 + pnorm(ur, nmean, nsd) - pnorm(ul, nmean, nsd) + pgpd(ur + epsilon, ur, sigmaur, xir))
  whichul = whichnonmiss[which(p[whichnonmiss] <= plul[whichnonmiss])]
  whichil = whichnonmiss[which((p[whichnonmiss] > plul[whichnonmiss]) & (p[whichnonmiss] < plbl[whichnonmiss]))]
  whichb = whichnonmiss[which((p[whichnonmiss] >= plbl[whichnonmiss]) & (p[whichnonmiss] <= plbr[whichnonmiss]))]
  whichir = whichnonmiss[which((p[whichnonmiss] > plbr[whichnonmiss]) & (p[whichnonmiss] < plur[whichnonmiss]))]
  whichur = whichnonmiss[which(p[whichnonmiss] >= plur[whichnonmiss])]
  nul = length(whichul)
  nil = length(whichil)
  nb = length(whichb)
  nir = length(whichir)
  nur = length(whichur)
  
  if (nul > 0) q[whichul] = -qgpd(1 - p[whichul]/kappa[whichul], -ul[whichul], sigmaul[whichul], xil[whichul])
  if (nb > 0) q[whichb] = qnorm(p[whichb]/kappa[whichb] - 1 + pnorm(ul[whichb], nmean[whichb], nsd[whichb]), nmean[whichb], nsd[whichb])
  if (nur > 0) {
    q[whichur] = qgpd(p[whichur]/kappa[whichur] - 1 - 
                      pnorm(ur[whichur], nmean[whichur], nsd[whichur]) + 
                      pnorm(ul[whichur], nmean[whichur], nsd[whichur]), ur[whichur], sigmaur[whichur], xir[whichur])
  }

  if (nil > 0) {
    for (i in 1:length(whichil)) {
      q[whichil[i]] = findqdmm(p[whichil[i]], 
                    ul[whichil[i]] + epsilon[whichil[i]] * c(-1, 1),
                    nmean[whichil[i]], nsd[whichil[i]], epsilon[whichil[i]], 
                    ul[whichil[i]], sigmaul[whichil[i]], xil[whichil[i]],
                    ur[whichil[i]], sigmaur[whichil[i]], xir[whichil[i]])
    }
  }
  if (nir > 0) {
    for (i in 1:length(whichir)) {
      q[whichir[i]] = findqdmm(p[whichir[i]],
                    ur[whichir[i]] + epsilon[whichir[i]] * c(-1, 1),
                    nmean[whichir[i]], nsd[whichir[i]], epsilon[whichir[i]], 
                    ul[whichir[i]], sigmaul[whichir[i]], xil[whichir[i]],
                    ur[whichir[i]], sigmaur[whichir[i]], xir[whichir[i]])
    }
  }
                  
  q
}

#' @export
#' @aliases itmgng ditmgng pitmgng qitmgng ritmgng
#' @rdname  itmgng

# random number generation for normal bulk with GPD's for upper and lower tails
ritmgng <- function(n = 1, nmean = 0, nsd = 1, epsilon = sd,
  ul = qnorm(0.1, nmean, nsd), sigmaul = nsd, xil = 0,
  ur = qnorm(0.9, nmean, nsd), sigmaur = nsd, xir = 0) {
  
  # Check properties of inputs
  check.n(n)
  check.param(nmean, allowvec = TRUE)
  check.posparam(nsd, allowvec = TRUE)
  check.posparam(epsilon, allowvec = TRUE, allowzero = TRUE)
  check.param(ul, allowvec = TRUE)
  check.posparam(sigmaul, allowvec = TRUE)
  check.param(xil, allowvec = TRUE)
  check.param(ur, allowvec = TRUE)
  check.posparam(sigmaur, allowvec = TRUE)
  check.param(xir, allowvec = TRUE)
  
  check.inputn(c(length(nmean), length(nsd), length(epsilon),
    length(ul), length(sigmaul), length(xil), length(ur), length(sigmaur), length(xir)))

  if (any(xil == 1) | any(xir == 1)) stop("shape cannot be 1")

  qitmgng(runif(n), nmean, nsd, epsilon, ul, sigmaul, xil, ur, sigmaur, xir)
}
