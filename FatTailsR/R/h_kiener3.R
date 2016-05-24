

#' @include g_kiener2.R



#' @title Asymmetric Kiener Distribution K3
#'
#' @description
#' Density, distribution function, quantile function, random generation
#' and additional formulae for asymmetric Kiener distribution K3.
#'
#' @param    x    vector of quantiles.
#' @param    q    vector of quantiles.
#' @param    m    numeric. The median.
#' @param    g    numeric. The scale parameter, preferably strictly positive.
#' @param    k    numeric. The tail parameter, preferably strictly positive.
#' @param    d    numeric. The distortion parameter between left and right tails.
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
#' Kiener distributions \code{K3(m, g, k, d, ...)} are distributions 
#' with asymmetrical left and right fat tails described by a global tail 
#' parameter \code{k} and a distortion parameter \code{d}. 
#' 
#' Distributions K3 (\code{\link{kiener3}}) 
#' with parameters \code{k} (kappa) and \code{d} (delta) and
#' distributions K4 (\code{\link{kiener4}})
#' with parameters \code{k} (kappa) and \code{e} (epsilon))
#' have been created to disantangle the parameters 
#' \code{a} (alpha) and \code{w} (omega) of distributions of 
#' distribution K2 (\code{\link{kiener2}}). 
#' The tiny difference between distributions K3 and K4 (\eqn{d = e/k}) 
#' has not yet been fully evaluated. Both should be tested at that moment.
#' 
#' \code{k} is the harmonic mean of \code{a} and \code{w} and represents a 
#' global tail parameter.
#'
#' \code{d} is a distortion parameter between the left tail parameter
#' \code{a} and the right tail parameter \code{w}.
#' It verifies the inequality: \eqn{-k < d < k} 
#' (whereas \code{e} of distribution K4 verifies \eqn{-1 < e < 1}).
#' The conversion functions (see \code{\link{aw2k}}) are:
#'
#' \deqn{1/k = (1/a + 1/w)/2 }
#' \deqn{  d = (-1/a + 1/w)/2 } 
#' \deqn{1/a = 1/k - d } 
#' \deqn{1/w = 1/k + d}
#'
#' \code{d} (and \code{e}) should be of the same sign than the skewness. 
#' A negative value \eqn{ d < 0 } implies \eqn{ a < w } and indicates a left  
#' tail heavier than the right tail. A positive value \eqn{ d > 0 } implies 
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
#' smaller than 1 or for very large absolute values of \code{d}. 
#' Hopefully, this case seldom happens in finance.
#'
#' \code{qkiener3} function is defined for p in (0, 1) by: 
#'   \deqn{ qkiener3(p, m, g, k, d) = 
#'               m + 2 * g * k * sinh(logit(p) / k) * exp(d * logit(p)) }
#'
#' \code{rkiener3} generates \code{n} random quantiles.
#'
#' In addition to the classical d, p, q, r functions, the prefixes 
#' dp, dq, l, dl, ql are also provided.
#'
#' \code{dpkiener3} is the density function calculated from the probability p. 
#' The formula is adapted from distribution K2. It is defined for p in (0, 1) by: 
#'   \deqn{ dpkiener3(p, m, g, k, d) = 
#'          p * (1 - p) / k / g / ( exp(-logit(p)/a)/a + exp(logit(p)/w)/w }
#' with \code{a} and \code{w} defined from \code{k} and \code{d} 
#' with the formula presented above.
#'
#' \code{dqkiener3} is the derivate of the quantile function calculated from 
#' the probability p. The formula is adapted from distribution K2. 
#' It is defined for p in (0, 1) by: 
#'   \deqn{ dqkiener3(p, m, g, k, d) = 
#'          k * g / p / (1 - p) * ( exp(-logit(p)/a)/a + exp(logit(p)/w)/w ) } 
#' with \code{a} and \code{w} defined above. 
#'
#' \code{lkiener3} function is estimated from a reverse optimization and can 
#' be (very) slow depending the number of points to estimate. Initialization 
#' is done with a symmetric distribution \code{\link{lkiener1}} 
#' of parameter \code{k} (thus \eqn{ d = 0}). Then optimization is performed  
#' to take into account the true value of \code{d}. 
#' The results can then be compared to the empirical probability logit(p).
#' WARNING: Results may become inconsistent when \code{k} is
#' smaller than 1 or for very large absolute values of \code{d}. 
#' Hopefully, this case seldom happens in finance.
#'
#' \code{dlkiener3} is the density function calculated from the logit of the 
#' probability lp = logit(p). The formula is adapted from distribution K2. 
#' it is defined for lp in (-Inf, +Inf) by: 
#'    \deqn{ dlkiener3(lp, m, g, k, d) = 
#'           p * (1 - p) / k / g / ( exp(-lp/a)/a + exp(lp/w)/w ) }
#' with \code{a} and \code{w} defined above. 
#'
#' \code{qlkiener3} is the quantile function calculated from the logit of the 
#' probability. It is defined for lp in (-Inf, +Inf) by: 
#'    \deqn{ qlkiener3(lp, m, g, k, d) = 
#'           m + 2 * g * k * sinh(lp / k) * exp(d * lp) }
#' 
#' \code{varkiener3} designates the Value a-risk and turns negative numbers 
#' into positive numbers with the following rule:
#'    \deqn{ varkiener3 <- if(p <= 0.5) { - qkiener3 } else { qkiener3 } }
#' Usual values in finance are \code{p = 0.01}, \code{p = 0.05}, \code{p = 0.95} and 
#' \code{p = 0.99}. \code{lower.tail = FALSE} uses \code{1-p} rather than {p}.
#' 
#' \code{ltmkiener3}, \code{rtmkiener3} and \code{eskiener3} are respectively the 
#' left tail mean, the right tail mean and the expected shortfall of the distribution 
#' (sometimes called average VaR, conditional VaR or tail VaR). 
#' Left tail mean is the integrale from \code{-Inf} to \code{p} of the quantile function 
#' \code{qkiener3} divided by \code{p}.
#' Right tail mean is the integrale from \code{p} to \code{+Inf} of the quantile function 
#' \code{qkiener3} divided by 1-p.
#' Expected shortfall turns negative numbers into positive numbers with the following rule:
#'    \deqn{ eskiener3 <- if(p <= 0.5) { - ltmkiener3 } else { rtmkiener3 } }
#' Usual values in finance are \code{p = 0.01}, \code{p = 0.025}, \code{p = 0.975} and 
#' \code{p = 0.99}. \code{lower.tail = FALSE} uses \code{1-p} rather than {p}.
#'
#' \code{dtmqkiener3} is the difference between the left tail mean and the quantile 
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
#' asymmetric Kiener distributions K2 and K4  
#' \code{\link{kiener2}}, \code{\link{kiener4}}, 
#' conversion functions \code{\link{aw2k}},  
#' estimation function \code{\link{fitkienerX}}, 
#' regression function \code{\link{regkienerLX}}.
#'
#' @examples
#' 
#' require(graphics)
#' 
#' ### Example 1
#' pp <- c(ppoints(11, a = 1), NA, NaN) ; pp
#' lp <- logit(pp) ; lp
#' qkiener3(  p = pp, m = 2, g = 1.5, k = aw2k(4, 6), d = aw2d(4, 6))
#' qlkiener3(lp = lp, m = 2, g = 1.5, k = aw2k(4, 6), d = aw2d(4, 6))
#' dpkiener3( p = pp, m = 2, g = 1.5, k = aw2k(4, 6), d = aw2d(4, 6))
#' dlkiener3(lp = lp, m = 2, g = 1.5, k = aw2k(4, 6), d = aw2d(4, 6))
#' dqkiener3( p = pp, m = 2, g = 1.5, k = aw2k(4, 6), d = aw2d(4, 6))
#' 
#' 
#' ### Example 2
#' k       <- 4.8
#' d       <- 0.042
#' set.seed(2014)
#' mainTC  <- paste("qkiener3(p, m = 0, g = 1, k = ", k, ", d = ", d, ")")
#' mainsum <- paste("cumulated qkiener3(p, m = 0, g = 1, k = ", k, ", d = ", d, ")")
#' T       <- 500
#' C       <- 4
#' TC      <- qkiener3(p = runif(T*C), m = 0, g = 1, k = k, d = d)
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
#' d     <- c(-0.1, -0.03, -0.01, 0.01, 0.03, 0.1) ; names(d) <- d
#' olty  <- c(2, 1, 2, 1, 2, 1, 1)
#' olwd  <- c(1, 1, 2, 2, 3, 3, 2)
#' ocol  <- c(2, 2, 4, 4, 3, 3, 1)
#' lleg  <- c("logit(0.999) = 6.9", "logit(0.99)   = 4.6", "logit(0.95)   = 2.9", 
#'            "logit(0.50)   = 0", "logit(0.05)   = -2.9", "logit(0.01)   = -4.6", 
#'            "logit(0.001) = -6.9  ")
#' op    <- par(mfrow=c(2,2), mgp=c(1.5,0.8,0), mar=c(3,3,2,1))
#' 
#' plot(x, pkiener3(x, k = 3.2, d = 0), type = "l", lwd = 3, ylim = c(0, 1),
#'      xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
#'      main = "pkiener3(q, m, g, k=3.2, d=...)")
#' for (i in 1:length(d)) lines(x, pkiener3(x, k = 3.2, d = d[i]), 
#'        lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' legend("topleft", title = expression(delta), legend = c(d, "0"), 
#'        cex = 0.7, inset = 0.02, lty = olty, lwd = olwd, col = ocol )
#' 
#' plot(x, dkiener3(x, k = 3.2, d = 0), type = "l", lwd = 3, ylim = c(0, 0.14),
#'      xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
#'      main = "dkiener3(q, m, g, k=3.2, d=...)")
#' for (i in 1:length(d)) lines(x, dkiener3(x, k = 3.2, d = d[i]), 
#'        lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' legend("topright", title = expression(delta), legend = c(d, "0"), 
#'        cex = 0.7, inset = 0.02, lty = olty, lwd = olwd, col = ocol )
#' 
#' plot(x, lkiener3(x, k = 3.2, d = 0), type = "l", lwd =3, ylim = c(-7.5, 7.5), 
#'      yaxt="n", xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
#'      main = "logit(pkiener3(q, m, g, k=3.2, d=...))")
#' axis(2, las=1, at=c(-6.9, -4.6, -2.9, 0, 2.9, 4.6, 6.9) )
#' for (i in 1:length(d)) lines(x, lkiener3(x, k = 3.2, d = d[i]),  
#'        lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' legend("topleft", legend = lleg, cex = 0.7, inset = 0.02 )
#' legend("bottomright", title = expression(delta), legend = c(d, "0"), 
#'        cex = 0.7, inset = 0.02, lty = c(olty), lwd = c(olwd), col = c(ocol) )
#' 
#' plot(x, dkiener3(x, k = 3.2, d = 0, log = TRUE), type = "l", lwd = 3, 
#'      ylim = c(-8, -1.5), xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
#'      main = "log(dkiener3(q, m, g, k=2, d=...))")
#' for (i in 1:length(d)) lines(x, dkiener3(x, k = 3.2, d = d[i], log=TRUE), 
#'        lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' legend("bottom", title = expression(delta), legend = c(d, "0"), 
#'        cex = 0.7, inset = 0.02, lty = olty, lwd = olwd, col = ocol )
#' ### End example 3
#' 
#' 
#' ### Example 4 (four plots: quantile, derivate, density and quantiles from p)
#' p     <- ppoints(199, a=0)
#' d     <- c(-0.1, -0.03, -0.01, 0.01, 0.03, 0.1) ; names(d) <- d
#' op    <- par(mfrow=c(2,2), mgp=c(1.5,0.8,0), mar=c(3,3,2,1))
#' 
#' plot(p, qlogis(p, scale = 2), type = "l", lwd = 2, xlim = c(0, 1),
#'       ylim = c(-15, 15), xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
#'      main = "qkiener3(p, m, g, k=3.2, d=...)")
#' for (i in 1:length(d)) lines(p, qkiener3(p, k = 3.2, d = d[i]), 
#'           lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' legend("topleft", title = expression(delta), legend = c(d, "qlogis(x/2)"), 
#'           inset = 0.02, lty = olty, lwd = olwd, col = ocol, cex = 0.7 )
#' 
#' plot(p, 2/p/(1-p), type = "l", lwd = 2, xlim = c(0, 1), ylim = c(0, 100),
#'      xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
#'      main = "dqkiener3(p, m, g, k=3.2, d=...)")
#' for (i in 1:length(d)) lines(p, dqkiener3(p, k = 3.2, d = d[i]), 
#'           lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' legend("top", title = expression(delta), legend = c(d, "p*(1-p)/2"), 
#'           inset = 0.02, lty = olty, lwd = olwd, col = ocol, cex = 0.7 )
#' 
#' plot(qlogis(p, scale = 2), p*(1-p)/2, type = "l", lwd = 2, xlim = c(-15, 15), 
#'      ylim = c(0, 0.14), xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
#'      main = "qkiener3, dpkiener3(p, m, g, k=3.2, d=...)")
#' for (i in 1:length(d)) { 
#'      lines(qkiener3(p, k = 3.2, d = d[i]), dpkiener3(p, k = 3.2, d = d[i]),
#'            lty = olty[i], lwd = olwd[i], col = ocol[i] ) }
#' legend("topleft", title = expression(delta), legend = c(d, "p*(1-p)/2"), 
#'           inset = 0.02, lty = olty, lwd = olwd, col = ocol, cex = 0.7 )
#' 
#' plot(qlogis(p, scale = 2), p, type = "l", lwd = 2, xlim = c(-15, 15), 
#'      ylim = c(0, 1), xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
#'      main = "inverse axis qkiener3(p, m, g, k=3.2, d=...)")
#' for (i in 1:length(d)) lines(qkiener3(p, k = 3.2, d = d[i]), p,
#'           lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' legend("topleft", title = expression(delta), legend = c(d, "qlogis(x/2)"), 
#'           inset = 0.02, lty = olty, lwd = olwd, col = ocol, cex = 0.7 )
#' ### End example 4
#' 
#' 
#' ### Example 5 (q and VaR, ltm, rtm, and ES)
#' pp <- c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 
#'         0.10, 0.20, 0.35, 0.5, 0.65, 0.80, 0.90,
#'         0.95, 0.975, 0.99, 0.995, 0.9975, 0.999)
#' m <- -10 ; g <- 1 ; k <- 4 ; d <- 0.06 
#' a <- dk2a(d, k) ; w <- dk2w(d, k) ; e <- dk2e(d, k)
#' round(c(m = m, g = g, a = a, k = k, w = w, d = d, e = e), 2) 
#' plot(qkiener3(  pp, m=m, k=k, d=d), pp, type ="b")
#' round(cbind(p = pp, "1-p" = 1-pp,
#' 	q   =   qkiener3(pp, m, g, k, d),
#' 	ltm = ltmkiener3(pp, m, g, k, d),
#' 	rtm = rtmkiener3(pp, m, g, k, d),
#' 	ES  =  eskiener3(pp, m, g, k, d),
#' 	VaR = varkiener3(pp, m, g, k, d)), 4)
#' round(kmean(c(m, g, k, d), model = "K3"), 4) # limit value for ltm and rtm
#' round(cbind(p = pp, "1-p" = 1-pp,
#' 	q   =   qkiener3(pp, m, g, k, d, lower.tail = FALSE),
#' 	ltm = ltmkiener3(pp, m, g, k, d, lower.tail = FALSE),
#' 	rtm = rtmkiener3(pp, m, g, k, d, lower.tail = FALSE),
#' 	ES  =  eskiener3(pp, m, g, k, d, lower.tail = FALSE),
#' 	VaR = varkiener3(pp, m, g, k, d, lower.tail = FALSE)), 4)
#' ### End example 5
#' 
#' 
#' @name kiener3
NULL

#' @export 
#' @rdname kiener3
dkiener3 <- function(x, m = 0, g = 1, k = 3.2, d = 0, log = FALSE) {
	lp <-  lkiener3(x,  m, g, k, d)
	v  <- dlkiener3(lp, m, g, k, d)
	if(log) return(log(v)) else return(v)
} 

#' @export
#' @rdname kiener3
pkiener3 <- function(q, m = 0, g = 1, k = 3.2, d = 0, 
                     lower.tail = TRUE, log.p = FALSE) {
	lp <- lkiener3(x = q, m, g, k, d)
	if(lower.tail) v <- invlogit(lp) else v <- 1 - invlogit(lp)
	if(log.p) return(log(v)) else return(v)
} 

#' @export
#' @rdname kiener3
qkiener3 <- function(p, m = 0, g = 1, k = 3.2, d = 0, 
                     lower.tail = TRUE, log.p = FALSE) {
	if(log.p) p <- exp(p) else p <- p
	if(lower.tail) p <- p else p <- 1-p
	v <- m + 2 * g * k * sinh(logit(p) / k) * exp(d * logit(p))
	return(v)
} 

#' @export
#' @rdname kiener3
rkiener3 <- function(n, m = 0, g = 1, k = 3.2, d = 0) {
	p <- runif(n)
	v <- qkiener3(p, m, g, k, d)
	return(v)
} 

#' @export
#' @rdname kiener3
dpkiener3 <- function(p, m = 0, g = 1, k = 3.2, d = 0, log = FALSE) {
	a <- kd2a(k, d)
	w <- kd2w(k, d)
	v <- p * (1 - p) / k / g / ( exp(-logit(p)/a)/a + exp(logit(p)/w)/w )
	if(log) return(log(v)) else return(v)
} 

#' @export
#' @rdname kiener3
dqkiener3 <- function(p, m = 0, g = 1, k = 3.2, d = 0, log = FALSE) {
# Compute dX/dp
	a <- kd2a(k, d)
	w <- kd2w(k, d)
	v <- k * g / p / (1 - p) * ( exp(-logit(p)/a)/a + exp(logit(p)/w)/w )
	if(log) return(log(v)) else return(v)
} 

#' @export
#' @rdname kiener3
lkiener3 <- function(x, m = 0, g = 1, k = 3.2, d = 0) { 
	lp.ini <- lkiener1(x, m, g, k)
	f      <- function(lp) sum( ( x - qlkiener3(lp, m, g, k, d) )^2 )
	lp.fin <- nlm(f, lp.ini)
	v      <- lp.fin$estimate
	return(v)
} 

#' @export
#' @rdname kiener3
dlkiener3 <- function(lp, m = 0, g = 1, k = 3.2, d = 0, log = FALSE) {
	p <- invlogit(lp)
	a <- kd2a(k, d)
	w <- kd2w(k, d)
	v <- p * (1 - p) / k / g / ( exp(-lp/a)/a + exp(lp/w)/w )
	if(log) return(log(v)) else return(v)
} # OK

#' @export
#' @rdname kiener3
qlkiener3 <- function(lp, m = 0, g = 1, k = 3.2, d = 0, lower.tail = TRUE ) {
	if(lower.tail) lp <- lp else lp <- -lp
	v  <- m + 2 * g * k * sinh(lp / k) * exp(d * lp)
	return(v)
} 

#' @export
#' @rdname kiener3
varkiener3 <- function(p, m = 0, g = 1, k = 3.2, d = 0, 
                      lower.tail = TRUE, log.p = FALSE) {
	p   <- if(log.p) {exp(p)} else {p}
	p   <- if(lower.tail) {p} else {1-p}
	va  <- p
	for (i in seq_along(p)) {
		va[i] <- ifelse(p[i] <= 0.5, 
					- qkiener3(p[i], m, g, k, d),
					  qkiener3(p[i], m, g, k, d))
	}
return(va)
}

#' @export
#' @rdname kiener3
ltmkiener3 <- function(p, m = 0, g = 1, k = 3.2, d = 0, 
                     lower.tail = TRUE, log.p = FALSE) {
	p  <- if(log.p) {exp(p)} else {p}
	ltm <- if (lower.tail) {
		m+g*k/p*(
			-pbeta(p,1+d-1/k, 1-d+1/k)*beta(1+d-1/k, 1-d+1/k)
			+pbeta(p,1+d+1/k, 1-d-1/k)*beta(1+d+1/k, 1-d-1/k))	
	} else {
		m+g*k/p*(
			-pbeta(p,1-d+1/k, 1+d-1/k)*beta(1-d+1/k, 1+d-1/k)
			+pbeta(p,1-d-1/k, 1+d+1/k)*beta(1-d-1/k, 1+d+1/k))
	}
return(ltm)
} 

#' @export
#' @rdname kiener3
rtmkiener3 <- function(p, m = 0, g = 1, k = 3.2, d = 0, 
                       lower.tail = TRUE, log.p = FALSE) {
	p   <- if(log.p) {exp(p)} else {p}
	rtm <- if (!lower.tail) {
		m+g*k/(1-p)*(
			-pbeta(1-p, 1+d-1/k, 1-d+1/k)*beta(1+d-1/k, 1-d+1/k)
			+pbeta(1-p, 1+d+1/k, 1-d-1/k)*beta(1+d+1/k, 1-d-1/k))	
	} else {
		m+g*k/(1-p)*(
			-pbeta(1-p, 1-d+1/k, 1+d-1/k)*beta(1-d+1/k, 1+d-1/k)
			+pbeta(1-p, 1-d-1/k, 1+d+1/k)*beta(1-d-1/k, 1+d+1/k))
	}
return(rtm)
}

#' @export
#' @rdname kiener3
dtmqkiener3 <- function(p, m = 0, g = 1, k = 3.2, d = 0, 
                        lower.tail = TRUE, log.p = FALSE) {
	dtmq <- p
	for (i in seq_along(p)) {
		dtmq[i] <- ifelse(p[i] <= 0.5, 
			ltmkiener3(p[i], m, g, k, d, lower.tail, log.p) 
			- qkiener3(p[i], m, g, k, d, lower.tail, log.p),
			rtmkiener3(p[i], m, g, k, d, lower.tail, log.p) 
			- qkiener3(p[i], m, g, k, d, lower.tail, log.p)) 
	}
return(dtmq)
}

#' @export
#' @rdname kiener3
eskiener3 <- function(p, m = 0, g = 1, k = 3.2, d = 0, 
                      lower.tail = TRUE, log.p = FALSE) {
	p   <- if(log.p) {exp(p)} else {p}
	p   <- if(lower.tail) {p} else {1-p}
	es  <- p
	for (i in seq_along(p)) {
		es[i] <- ifelse(p[i] <= 0.5, 
					- ltmkiener3(p[i], m, g, k, d),
					  rtmkiener3(p[i], m, g, k, d))
	}
return(es)
}


