

#' @include h_kiener3.R



#' @title Asymmetric Kiener Distribution K4
#'
#' @description
#' Density, distribution function, quantile function, random generation
#' and additional formulae for asymmetric Kiener distribution K4.
#'
#' @param    x    vector of quantiles.
#' @param    q    vector of quantiles.
#' @param    m    numeric. The median.
#' @param    g    numeric. The scale parameter, preferably strictly positive.
#' @param    k    numeric. The tail parameter, preferably strictly positive.
#' @param    e    numeric. The eccentricity parameter between left and right tails.
#' @param    p    vector of probabilities.
#' @param    lp   vector of logit of probabilities.
#' @param    n    number of observations. If length(n) > 1, the length is  
#'                taken to be the number required.
#' @param    log           logical. If TRUE, densities are given in log scale.
#' @param    lower.tail    logical. If TRUE, use p. If FALSE, use 1-p.
#' @param    log.p         logical. If TRUE, probabilities p are given as log(p).
#'
#' @details
#' Kiener distributions use the following parameters, some of them being redundant. 
#' See \code{\link{aw2k}} and \code{\link{pk2pk}} for the formulas and 
#' the conversion between parameters:
#' \itemize{
#'   \item{ \code{m} (mu) is the median of the distribution,. }
#'   \item{ \code{g} (gamma) is the scale parameter. }
#'   \item{ \code{a} (alpha) is the left tail parameter. } 
#'   \item{ \code{k} (kappa) is the harmonic mean of \code{a} and \code{w} 
#'          and describes a global tail parameter. }
#'   \item{ \code{w} (omega) is the right tail parameter. } 
#'   \item{ \code{d} (delta) is the distortion parameter. }
#'   \item{ \code{e} (epsilon) is the eccentricity parameter. }
#' }
#' 
#' Kiener distributions \code{K4(m, g, k, e, ...)} are distributions 
#' with asymmetrical left and right fat tails described by a global tail 
#' parameter \code{k} and an eccentricity parameter \code{e}. 
#'  
#' Distributions K3 (\code{\link{kiener3}}) 
#' with parameters \code{k} (kappa) and \code{d} (delta) and
#' distributions K4 (\code{\link{kiener4}})
#' with parameters \code{k} (kappa) and \code{e} (epsilon))
#' have been created to disantangle the parameters 
#' \code{a} (alpha) and \code{w} (omega) of distributions K2
#' (\code{\link{kiener2}}). 
#' The tiny difference between distributions K3 and K4 (\eqn{d = e/k}) 
#' has not yet been fully evaluated. Both should be tested at that moment.
#' 
#' \code{k} is the harmonic mean of \code{a} and \code{w} and represents a 
#' global tail parameter.
#'
#' \code{e} is an eccentricity parameter between the left tail parameter
#' \code{a} and the right tail parameter \code{w}.
#' It verifies the inequality: \eqn{-1 < e < 1} 
#' (whereas \code{d} of distribution K3 verifies \eqn{-k < d < k}).
#' The conversion functions (see \code{\link{aw2k}}) are:
#'
#' \deqn{1/k = (1/a + 1/w)/2 }
#' \deqn{  e = (a - w)/(a + w) }
#' \deqn{  a = k/(1 - e) }
#' \deqn{  w = k/(1 + e) }
#'
#' \code{e} (and \code{d}) should be of the same sign than the skewness. 
#' A negative value \eqn{ e < 0 } implies \eqn{ a < w } and indicates a left 
#' tail heavier than the right tail. A positive value \eqn{ e > 0 } implies 
#' \eqn{ a > w } and a right tail heavier than the left tail.  
#'
#' \code{m} is the median of the distribution. \code{g} is the scale parameter 
#' and the inverse of the density at the median: \eqn{ g = 1 / 8 / f(m) }.
#' As a first estimate, it is approximatively one fourth of the standard 
#' deviation \eqn{ g  \approx \sigma / 4 } but is independant from it.
#' 
#' The d, p functions have no explicit forms. They are provided here for 
#' convenience. They are estimated from a reverse optimization on the quantile 
#' function and can be (very) slow, depending the number of points to estimate. 
#' We recommand to use the quantile function as far as possible. 
#' WARNING: Results may become inconsistent when \code{k} is
#' smaller than 1 or for very large absolute values of \code{e}. 
#' Hopefully, these cases seldom happen in finance.
#'
#' \code{qkiener4} function is defined for p in (0, 1) by: 
#'   \deqn{ qkiener4(p, m, g, k, e) = 
#'               m + 2 * g * k * sinh(logit(p) / k) * exp(e / k * logit(p)) }
#'
#' \code{rkiener4} generates \code{n} random quantiles.
#'
#' In addition to the classical d, p, q, r functions, the prefixes 
#' dp, dq, l, dl, ql are also provided.
#'
#' \code{dpkiener4} is the density function calculated from the probability p. 
#' The formula is adapted from distribution K2. It is defined for p in (0, 1) by: 
#'   \deqn{ dpkiener4(p, m, g, k, e) = 
#'          p * (1 - p) / k / g / ( exp(-logit(p)/a)/a + exp(logit(p)/w)/w }
#' with \code{a} and \code{w} defined from \code{k} and \code{e}.
#'
#' \code{dqkiener4} is the derivate of the quantile function calculated from 
#' the probability p. The formula is adapted from distribution K2. 
#' It is defined for p in (0, 1) by: 
#'   \deqn{ dqkiener4(p, m, g, k, e) = 
#'          k * g / p / (1 - p) * ( exp(-logit(p)/a)/a + exp(logit(p)/w)/w ) } 
#' with \code{a} and \code{w} defined with the formula presented above. 
#'
#' \code{lkiener4} function is estimated from a reverse optimization and can 
#' be (very) slow depending the number of points to estimate. Initialization 
#' is done with a symmetric distribution \code{\link{lkiener1}} 
#' of parameter \code{k} (thus \eqn{ e = 0}). Then optimization is performed  
#' to take into account the true value of \code{e}.  
#' The results can then be compared to the empirical probability logit(p).
#' WARNING: Results may become inconsistent when \code{k} is
#' smaller than 1 or for very large absolute values of \code{e}. 
#' Hopefully, these cases seldom happen in finance.
#'
#' \code{dlkiener4} is the density function calculated from the logit of the 
#' probability lp = logit(p). The formula is adapted from distribution K2.
#' it is defined for lp in (-Inf, +Inf) by: 
#'    \deqn{ dlkiener4(lp, m, g, k, e) = 
#'           p * (1 - p) / k / g / ( exp(-lp/a)/a + exp(lp/w)/w ) }
#' with \code{a} and \code{w} defined above. 
#'
#' \code{qlkiener4} is the quantile function calculated from the logit of the 
#' probability. It is defined for lp in (-Inf, +Inf) by: 
#'    \deqn{ qlkiener4(lp, m, g, k, e) = 
#'           m + 2 * g * k * sinh(lp / k) * exp(e / k * lp) }
#' 
#' \code{varkiener4} designates the Value a-risk and turns negative numbers 
#' into positive numbers with the following rule:
#'    \deqn{ varkiener4 <- if(p <= 0.5) { - qkiener4 } else { qkiener4 } }
#' Usual values in finance are \code{p = 0.01}, \code{p = 0.05}, \code{p = 0.95} and 
#' \code{p = 0.99}. \code{lower.tail = FALSE} uses \code{1-p} rather than {p}.
#' 
#' \code{ltmkiener4}, \code{rtmkiener4} and \code{eskiener4} are respectively the 
#' left tail mean, the right tail mean and the expected shortfall of the distribution 
#' (sometimes called average VaR, conditional VaR or tail VaR). 
#' Left tail mean is the integrale from \code{-Inf} to \code{p} of the quantile function 
#' \code{qkiener4} divided by \code{p}.
#' Right tail mean is the integrale from \code{p} to \code{+Inf} of the quantile function 
#' \code{qkiener4} divided by 1-p.
#' Expected shortfall turns negative numbers into positive numbers with the following rule:
#'    \deqn{ eskiener4 <- if(p <= 0.5) { - ltmkiener4 } else { rtmkiener4 } }
#' Usual values in finance are \code{p = 0.01}, \code{p = 0.025}, \code{p = 0.975} and 
#' \code{p = 0.99}. \code{lower.tail = FALSE} uses \code{1-p} rather than {p}.
#'
#' \code{dtmqkiener4} is the difference between the left tail mean and the quantile 
#' when (p <= 0.5) and the difference between the right tail mean and the quantile 
#' when (p > 0.5). It is in quantile unit and is an indirect measure of the tail curvature.
#' 
#' @references
#' P. Kiener, Explicit models for bilateral fat-tailed distributions and 
#' applications in finance with the package FatTailsR, 8th R/Rmetrics Workshop 
#' and Summer School, Paris, 27 June 2014. Download it from:  
#' \url{http://www.inmodelia.com/exemples/2014-0627-Rmetrics-Kiener-en.pdf}
#'
#' P. Kiener, Fat tail analysis and package FatTailsR, 
#' 9th R/Rmetrics Workshop and Summer School, Zurich, 27 June 2015. 
#' Download it from: 
#' \url{http://www.inmodelia.com/exemples/2015-0627-Rmetrics-Kiener-en.pdf}
#' 
#' C. Acerbi, D. Tasche, Expected shortfall: a natural coherent alternative to 
#' Value at Risk, 9 May 2001. Download it from: 
#' \url{http://www.bis.org/bcbs/ca/acertasc.pdf}
#' 
#' @seealso 
#' Symmetric Kiener distribution K1 \code{\link{kiener1}}, 
#' asymmetric Kiener distributions K2 and K3 
#' \code{\link{kiener2}}, \code{\link{kiener3}}, 
#' conversion functions \code{\link{aw2k}}, 
#' estimation function \code{\link{fitkienerX}},
#  regression function \code{\link{regkienerLX}}.
#'
#' @examples
#' 
#' require(graphics)
#' 
#' ### Example 1
#' pp <- c(ppoints(11, a = 1), NA, NaN) ; pp
#' lp <- logit(pp) ; lp
#' qkiener4(  p = pp, m = 2, g = 1.5, k = aw2k(4, 6), e = aw2e(4, 6))
#' qlkiener4(lp = lp, m = 2, g = 1.5, k = aw2k(4, 6), e = aw2e(4, 6))
#' dpkiener4( p = pp, m = 2, g = 1.5, k = aw2k(4, 6), e = aw2e(4, 6))
#' dlkiener4(lp = lp, m = 2, g = 1.5, k = aw2k(4, 6), e = aw2e(4, 6))
#' dqkiener4( p = pp, m = 2, g = 1.5, k = aw2k(4, 6), e = aw2e(4, 6))
#' 
#' 
#' ### Example 2
#' k       <- 4.8
#' e       <- 0.2
#' set.seed(2014)
#' mainTC  <- paste("qkiener4(p, m = 0, g = 1, k = ", k, ", e = ", e, ")")
#' mainsum <- paste("cumulated qkiener4(p, m = 0, g = 1, k = ", k, ", e = ", e, ")")
#' T       <- 500
#' C       <- 4
#' TC      <- qkiener4(p = runif(T*C), m = 0, g = 1, k = k, e = e)
#' matTC   <- matrix(TC, nrow = T, ncol = C, dimnames = list(1:T, letters[1:C]))
#' head(matTC)
#' plot.ts(matTC, main = mainTC)
#' #
#' matsum  <- apply(matTC, MARGIN=2, cumsum) 
#' head(matsum)
#' plot.ts(matsum, plot.type = "single", main = mainsum)
#' ### End example 2
#' 
#' 
#' ### Example 3 (four plots: probability, density, logit, logdensity)
#' x     <- q  <- seq(-15, 15, length.out=101)
#' k     <- 3.2
#' e     <- c(-0.3, -0.15, -0.07, 0.07, 0.15, 0.30) ; names(e) <- e
#' olty  <- c(2, 1, 2, 1, 2, 1, 1)
#' olwd  <- c(1, 1, 2, 2, 3, 3, 2)
#' ocol  <- c(2, 2, 4, 4, 3, 3, 1)
#' lleg  <- c("logit(0.999) = 6.9", "logit(0.99)   = 4.6", "logit(0.95)   = 2.9", 
#'            "logit(0.50)   = 0", "logit(0.05)   = -2.9", "logit(0.01)   = -4.6", 
#'            "logit(0.001) = -6.9  ")
#' op    <- par(mfrow=c(2,2), mgp=c(1.5,0.8,0), mar=c(3,3,2,1))
#' 
#' plot(x, pkiener4(x, k = 3.2, e = 0), type = "l", lwd = 3, ylim = c(0, 1),
#'      xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
#'      main = "pkiener4(q, m, g, k=3.2, e=...)")
#' for (i in 1:length(e)) lines(x, pkiener4(x, k = 3.2, e = e[i]), 
#'        lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' legend("topleft", title = expression(epsilon), legend = c(e, "0"), 
#'        cex = 0.7, inset = 0.02, lty = olty, lwd = olwd, col = ocol )
#' 
#' plot(x, dkiener4(x, k = 3.2, e = 0), type = "l", lwd = 3, ylim = c(0, 0.14),
#'      xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
#'      main = "dkiener4(q, m, g, k=3.2, e=...)")
#' for (i in 1:length(e)) lines(x, dkiener4(x, k = 3.2, e = e[i]), 
#'        lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' legend("topright", title = expression(epsilon), legend = c(e, "0"), 
#'        cex = 0.7, inset = 0.02, lty = olty, lwd = olwd, col = ocol )
#' 
#' plot(x, lkiener4(x, k = 3.2, e = 0), type = "l", lwd =3, ylim = c(-7.5, 7.5), 
#'      yaxt="n", xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
#'      main = "logit(pkiener4(q, m, g, k=3.2, e=...))")
#' axis(2, las=1, at=c(-6.9, -4.6, -2.9, 0, 2.9, 4.6, 6.9) )
#' for (i in 1:length(e)) lines(x, lkiener4(x, k = 3.2, e = e[i]),  
#'        lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' legend("topleft", legend = lleg, cex = 0.7, inset = 0.02 )
#' legend("bottomright", title = expression(epsilon), legend = c(e, "0"), 
#'        cex = 0.7, inset = 0.02, lty = c(olty), lwd = c(olwd), col = c(ocol) )
#' 
#' plot(x, dkiener4(x, k = 3.2, e = 0, log = TRUE), type = "l", lwd = 3, 
#'      ylim = c(-8, -1.5), xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
#'      main = "log(dkiener4(q, m, g, k=2, e=...))")
#' for (i in 1:length(e)) lines(x, dkiener4(x, k = 3.2, e = e[i], log=TRUE), 
#'        lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' legend("bottom", title = expression(epsilon), legend = c(e, "0"), 
#'        cex = 0.7, inset = 0.02, lty = olty, lwd = olwd, col = ocol )
#' ### End example 3
#' 
#' 
#' ### Example 4 (four plots: quantile, derivate, density and quantiles from p)
#' p     <- ppoints(199, a=0)
#' e     <- c(-0.3, -0.15, -0.07, 0.07, 0.15, 0.30) ; names(e) <- e
#' op    <- par(mfrow=c(2,2), mgp=c(1.5,0.8,0), mar=c(3,3,2,1))
#' 
#' plot(p, qlogis(p, scale = 2), type = "l", lwd = 2, xlim = c(0, 1),
#'      ylim = c(-15, 15), xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
#'      main = "qkiener4(p, m, g, k=3.2, e=...)")
#' for (i in 1:length(e)) lines(p, qkiener4(p, k = 3.2, e = e[i]), 
#'           lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' legend("topleft", title = expression(epsilon), legend = c(e, "qlogis(x/2)"), 
#'         inset = 0.02, lty = olty, lwd = olwd, col = ocol, cex = 0.7 )
#' 
#' plot(p, 2/p/(1-p), type = "l", lwd = 2, xlim = c(0, 1), ylim = c(0, 100),
#'      xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
#'      main = "dqkiener4(p, m, g, k=3.2, e=...)")
#' for (i in 1:length(e)) lines(p, dqkiener4(p, k = 3.2, e = e[i]), 
#'           lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' legend("top", title = expression(epsilon), legend = c(e, "p*(1-p)/2"), 
#'           inset = 0.02, lty = olty, lwd = olwd, col = ocol, cex = 0.7 )
#' 
#' plot(qlogis(p, scale = 2), p*(1-p)/2, type = "l", lwd = 2, xlim = c(-15, 15), 
#'      ylim = c(0, 0.14), xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
#'      main = "qkiener4, dpkiener4(p, m, g, k=3.2, e=...)")
#' for (i in 1:length(e)) { 
#'      lines(qkiener4(p, k = 3.2, e = e[i]), dpkiener4(p, k = 3.2, e = e[i]),
#'            lty = olty[i], lwd = olwd[i], col = ocol[i] ) }
#' legend("topleft", title = expression(epsilon), legend = c(e, "p*(1-p)/2"), 
#'           inset = 0.02, lty = olty, lwd = olwd, col = ocol, cex = 0.7 )
#' 
#' plot(qlogis(p, scale = 2), p, type = "l", lwd = 2, xlim = c(-15, 15), 
#'      ylim = c(0, 1), xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
#'      main = "inverse axis qkiener4(p, m, g, k=3.2, e=...)")
#' for (i in 1:length(e)) lines(qkiener4(p, k = 3.2, e = e[i]), p,
#'           lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' legend("topleft", title = expression(epsilon), legend = c(e, "qlogis(x/2)"), 
#'           inset = 0.02, lty = olty, lwd = olwd, col = ocol, cex = 0.7 )
#' ### End example 4
#' 
#' 
### Example 5 (q and VaR, ltm, rtm, and ES)
#' pp <- c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 
#'         0.10, 0.20, 0.35, 0.5, 0.65, 0.80, 0.90,
#'         0.95, 0.975, 0.99, 0.995, 0.9975, 0.999)
#' m <- -5 ; g <- 1 ; k <- 4 ; e = -0.20
#' a <- ek2a(e, k) ; w <- ek2w(e, k) ; d <- ek2d(e, k) 
#' round(c(m = m, g = g, a = a, k = k, w = w, d = d, e = e), 2) 
#' plot(qkiener4(pp, m, g, k, e), pp, type = "b")
#' round(cbind(p = pp, "1-p" = 1-pp,
#' 	q   =   qkiener4(pp, m, g, k, e),
#' 	ltm = ltmkiener4(pp, m, g, k, e),
#' 	rtm = rtmkiener4(pp, m, g, k, e),
#' 	ES  =  eskiener4(pp, m, g, k, e),
#' 	VaR = varkiener4(pp, m, g, k, e)), 4)
#' round(kmean(c(m, g, k, e), model = "K4"), 4) # limit value for ltm and rtm
#' round(cbind(p = pp, "1-p" = 1-pp, 
#' 	q   =   qkiener4(pp, m, g, k, e, lower.tail = FALSE), 
#' 	ltm = ltmkiener4(pp, m, g, k, e, lower.tail = FALSE), 
#' 	rtm = rtmkiener4(pp, m, g, k, e, lower.tail = FALSE), 
#' 	ES  =  eskiener4(pp, m, g, k, e, lower.tail = FALSE), 
#' 	VaR = varkiener4(pp, m, g, k, e, lower.tail = FALSE)), 4)
#' ### End example 5
#' 
#' 
#' @name kiener4
NULL

#' @export 
#' @rdname kiener4
dkiener4 <- function(x, m = 0, g = 1, k = 3.2, e = 0, log = FALSE) {
	lp <-  lkiener4(x,  m, g, k, e)
	v  <- dlkiener4(lp, m, g, k, e)
	if(log) return(log(v)) else return(v)
}

#' @export
#' @rdname kiener4
pkiener4 <- function(q, m = 0, g = 1, k = 3.2, e = 0, 
                     lower.tail = TRUE, log.p = FALSE) {
	lp <- lkiener4(x = q, m, g, k, e)
	if(lower.tail) v <- invlogit(lp) else v <- 1 - invlogit(lp)
	if(log.p) return(log(v)) else return(v)
}

#' @export
#' @rdname kiener4
qkiener4 <- function(p, m = 0, g = 1, k = 3.2, e = 0, 
                     lower.tail = TRUE, log.p = FALSE) {
	if(log.p) p <- exp(p) else p <- p
	if(lower.tail) p <- p else p <- 1-p
	v <- m + 2 * g * k * sinh(logit(p) / k) * exp(e / k * logit(p))
	return(v)
}

#' @export
#' @rdname kiener4
rkiener4 <- function(n, m = 0, g = 1, k = 3.2, e = 0) {
	p <- runif(n)
	v <- qkiener4(p, m, g, k, e)
	return(v)
}

#' @export
#' @rdname kiener4
dpkiener4 <- function(p, m = 0, g = 1, k = 3.2, e = 0, log = FALSE) {
	a <- ke2a(k, e)
	w <- ke2w(k, e)
	v <- p * (1 - p) / k / g / ( exp(-logit(p)/a)/a + exp(logit(p)/w)/w )
	if(log) return(log(v)) else return(v)
}

#' @export
#' @rdname kiener4
dqkiener4 <- function(p, m = 0, g = 1, k = 3.2, e = 0, log = FALSE) {
# Compute dX/dp
	a <- ke2a(k, e)
	w <- ke2w(k, e)
	v <- k * g / p / (1 - p) * ( exp(-logit(p)/a)/a + exp(logit(p)/w)/w )
	if(log) return(log(v)) else return(v)
}

## voir fonction NLSstClosestX
#' @export
#' @rdname kiener4
lkiener4 <- function(x, m = 0, g = 1, k = 3.2, e = 0) { 
	lp.ini <- lkiener1(x, m, g, k)
	f      <- function(lp) sum( ( x - qlkiener4(lp, m, g, k, e) )^2 )
	lp.fin <- nlm(f, lp.ini)
	v      <- lp.fin$estimate
	return(v)
} 

#' @export
#' @rdname kiener4
dlkiener4 <- function(lp, m = 0, g = 1, k = 3.2, e = 0, log = FALSE) {
	p <- invlogit(lp)
	a <- ke2a(k, e)
	w <- ke2w(k, e)
	v <- p * (1 - p) / k / g / ( exp(-lp/a)/a + exp(lp/w)/w )
	if(log) return(log(v)) else return(v)
} 

#' @export
#' @rdname kiener4
qlkiener4 <- function(lp, m = 0, g = 1, k = 3.2, e = 0, lower.tail = TRUE ) {
	if(lower.tail) lp <- lp else lp <- -lp
		v <- m + 2 * g * k * sinh(lp / k) * exp(e / k * lp)
	return(v)
} 

#' @export
#' @rdname kiener4
varkiener4 <- function(p, m = 0, g = 1, k = 3.2, e = 0, 
                      lower.tail = TRUE, log.p = FALSE) {
	p   <- if(log.p) {exp(p)} else {p}
	p   <- if(lower.tail) {p} else {1-p}
	va  <- p
	for (i in seq_along(p)) {
		va[i] <- ifelse(p[i] <= 0.5, 
					- qkiener4(p[i], m, g, k, e),
					  qkiener4(p[i], m, g, k, e))
	}
return(va)
}

#' @export
#' @rdname kiener4
ltmkiener4 <- function(p, m = 0, g = 1, k = 3.2, e = 0, 
                     lower.tail = TRUE, log.p = FALSE) {
	p <- if(log.p) {exp(p)} else {p}
	ltm <- if (lower.tail) {
		m+g*k/p*(
	        -pbeta(p, 1-(1-e)/k, 1+(1-e)/k) * beta(1-(1-e)/k, 1+(1-e)/k) 
	        +pbeta(p, 1+(1+e)/k, 1-(1+e)/k) * beta(1+(1+e)/k, 1-(1+e)/k))
	} else {
		m+g*k/p*(
			-pbeta(p, 1+(1-e)/k, 1-(1-e)/k)*beta(1+(1-e)/k, 1-(1-e)/k)
			+pbeta(p, 1-(1+e)/k, 1+(1+e)/k)*beta(1-(1+e)/k, 1+(1+e)/k))
	}
return(ltm)
}

#' @export
#' @rdname kiener4
rtmkiener4 <- function(p, m = 0, g = 1, k = 3.2, e = 0, 
                       lower.tail = TRUE, log.p = FALSE) {
	p   <- if(log.p) {exp(p)} else {p}
	rtm <- if (!lower.tail) {
		m+g*k/(1-p)*(
			-pbeta(1-p, 1-(1-e)/k, 1+(1-e)/k)*beta(1-(1-e)/k, 1+(1-e)/k)
			+pbeta(1-p, 1+(1+e)/k, 1-(1+e)/k)*beta(1+(1+e)/k, 1-(1+e)/k))	
	} else {
		m+g*k/(1-p)*(
			-pbeta(1-p, 1+(1-e)/k, 1-(1-e)/k)*beta(1+(1-e)/k, 1-(1-e)/k)
			+pbeta(1-p, 1-(1+e)/k, 1+(1+e)/k)*beta(1-(1+e)/k, 1+(1+e)/k))
	}
return(rtm)
}

#' @export
#' @rdname kiener4
dtmqkiener4 <- function(p, m = 0, g = 1, k = 3.2, e = 0, 
                      lower.tail = TRUE, log.p = FALSE) {
	dtmq <- p
	for (i in seq_along(p)) {
		dtmq[i] <- ifelse(p[i] <= 0.5, 
			ltmkiener4(p[i], m, g, k, e, lower.tail, log.p) 
			- qkiener4(p[i], m, g, k, e, lower.tail, log.p),
			rtmkiener4(p[i], m, g, k, e, lower.tail, log.p) 
			- qkiener4(p[i], m, g, k, e, lower.tail, log.p)) 
	}
return(dtmq)
}

#' @export
#' @rdname kiener4
eskiener4 <- function(p, m = 0, g = 1, k = 3.2, e = 0, 
                      lower.tail = TRUE, log.p = FALSE) {
	p   <- if(log.p) {exp(p)} else {p}
	p   <- if(lower.tail) {p} else {1-p}
	es  <- p
	for (i in seq_along(p)) {
		es[i] <- ifelse(p[i] <= 0.5, 
					- ltmkiener4(p[i], m, g, k, e),
					  rtmkiener4(p[i], m, g, k, e))
	}
return(es)
}

