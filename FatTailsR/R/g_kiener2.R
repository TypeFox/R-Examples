

#' @include f_kiener1.R



#' @title Asymmetric Kiener Distribution K2
#'
#' @description
#' Density, distribution function, quantile function, random generation
#' and additional formulae for asymmetric Kiener distribution K2.
#'
#' @param    x    vector of quantiles.
#' @param    q    vector of quantiles.
#' @param    m    numeric. The median.
#' @param    g    numeric. The scale parameter, preferably strictly positive.
#' @param    a    numeric. The left tail parameter, preferably strictly positive.
#' @param    w    numeric. The right tail parameter, preferably strictly positive.
#' @param    p    vector of probabilities.
#' @param    lp   vector of logit of probabilities.
#' @param    n    number of observations. If length(n) > 1, the length is  
#'                                  taken to be the number required.
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
#' Kiener distributions \code{K2(m, g, a, w)} are distributions 
#' with asymmetrical left 
#' and right fat tails described by the parameters \code{a} (alpha) for 
#' the left tail and \code{w} (omega) for the right tail. These parameters 
#' correspond to the power exponent that appear in Pareto formula and 
#' Karamata theorems. 
#'
#' As \code{a} and \code{w} are highly correlated, the use of Kiener distributions
#' (\code{K3(..., k, d)} K4 (\code{K4(..., k, e)} is an alternate solution.
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
#' WARNING: Results may become inconsistent when \code{a} or \code{w} are
#' smaller than 1. Hopefully, this case seldom happens in finance.
#'
#' \code{qkiener2} function is defined for p in (0, 1) by: 
#'   \deqn{ qkiener2(p, m, g, a, w) = 
#'                   m + g * k * (- exp(-logit(p)/a) + exp(logit(p)/w) ) }
#' where k is the harmonic mean of the tail parameters \code{a} and \code{w} 
#' calculated by \eqn{k = aw2k(a, w)}.
#'
#' \code{rkiener2} generates \code{n} random quantiles.
#'
#' In addition to the classical d, p, q, r functions, the prefixes 
#' dp, dq, l, dl, ql are also provided.
#'
#' \code{dpkiener2} is the density function calculated from the probability p. 
#' It is defined for p in (0, 1) by: 
#'   \deqn{ dpkiener2(p, m, g, a, w) = 
#'          p * (1 - p) / k / g / ( exp(-logit(p)/a)/a + exp(logit(p)/w)/w }
#'
#' \code{dqkiener2} is the derivate of the quantile function calculated from 
#' the probability p. It is defined for p in (0, 1) by: 
#'   \deqn{ dqkiener2(p, m, g, a, w) = 
#'          k * g / p / (1 - p) * ( exp(-logit(p)/a)/a + exp(logit(p)/w)/w ) }
#'
#' \code{lkiener2} function is estimated from a reverse optimization and can 
#' be (very) slow depending the number of points to estimate. Initialization 
#' is done by assuming a symmetric distribution \code{\link{lkiener1}} 
#' around the harmonic mean \code{k}, then optimization is performed to 
#' take into account the true values \code{a} and \code{w}. 
#' The result can be then compared to the empirical probability logit(p). 
#' WARNING: Results may become inconsistent when \code{a} or \code{w} are
#' smaller than 1. Hopefully, this case seldom happens in finance.
#'
#' \code{dlkiener2} is the density function calculated from the logit of the 
#' probability lp = logit(p).  
#' it is defined for lp in (-Inf, +Inf) by: 
#'    \deqn{ dlkiener2(lp, m, g, a, w) = 
#'           p * (1 - p) / k / g / ( exp(-lp/a)/a + exp(lp/w)/w ) }
#'
#' \code{qlkiener2} is the quantile function calculated from the logit of the 
#' probability. It is defined for lp in (-Inf, +Inf) by: 
#'    \deqn{ qlkiener2(lp, m, g, a, w) = 
#'           m + g * k * ( - exp(-lp/a) + exp(lp/w) ) }
#' 
#' \code{varkiener2} designates the Value a-risk and turns negative numbers 
#' into positive numbers with the following rule:
#'    \deqn{ varkiener2 <- if(p <= 0.5) { - qkiener2 } else { qkiener2 } }
#' Usual values in finance are \code{p = 0.01}, \code{p = 0.05}, \code{p = 0.95} and 
#' \code{p = 0.99}. \code{lower.tail = FALSE} uses \code{1-p} rather than {p}.
#' 
#' \code{ltmkiener2}, \code{rtmkiener2} and \code{eskiener2} are respectively the 
#' left tail mean, the right tail mean and the expected shortfall of the distribution 
#' (sometimes called average VaR, conditional VaR or tail VaR). 
#' Left tail mean is the integrale from \code{-Inf} to \code{p} of the quantile function 
#' \code{qkiener2} divided by \code{p}.
#' Right tail mean is the integrale from \code{p} to \code{+Inf} of the quantile function 
#' \code{qkiener2} divided by 1-p.
#' Expected shortfall turns negative numbers into positive numbers with the following rule:
#'    \deqn{ eskiener2 <- if(p <= 0.5) { - ltmkiener2 } else { rtmkiener2 } }
#' Usual values in finance are \code{p = 0.01}, \code{p = 0.025}, \code{p = 0.975} and 
#' \code{p = 0.99}. \code{lower.tail = FALSE} uses \code{1-p} rather than {p}.
#'
#' \code{dtmqkiener2} is the difference between the left tail mean and the quantile 
#' when (p <= 0.5) and the difference between the right tail mean and the quantile 
#' when (p > 0.5). It is in quantile unit and is an indirect measure of the tail curvature.
#' 
#' @references
#' P. Kiener, Explicit models for bilateral fat-tailed distributions and 
#' applications in finance with the package FatTailsR, 8th R/Rmetrics Workshop 
#' and Summer School, Paris, 27 June 2014.  Download it from:  
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
#' asymmetric Kiener distributions K3 and K4 
#' \code{\link{kiener3}}, \code{\link{kiener4}}, 
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
#' qkiener2(  p = pp, m = 2, g = 1.5, a = 4, w = 6)
#' qkiener2(  p = pp, m = 2, g = 1.5, a = 4, w = 6)
#' qlkiener2(lp = lp, m = 2, g = 1.5, a = 4, w = 6)
#' dpkiener2( p = pp, m = 2, g = 1.5, a = 4, w = 6)
#' dlkiener2(lp = lp, m = 2, g = 1.5, a = 4, w = 6)
#' dqkiener2( p = pp, m = 2, g = 1.5, a = 4, w = 6)
#' 
#' 
#' ### Example 2
#' a       <- 6
#' w       <- 4
#' set.seed(2014)
#' mainTC  <- paste("qkiener2(p, m = 0, g = 1, a = ", a, ", w = ", w, ")")
#' mainsum <- paste("cumulated qkiener2(p, m = 0, g = 1, a = ", a, ", w = ", w, ")")
#' T       <- 500
#' C       <- 4
#' TC      <- qkiener2(p = runif(T*C), m = 0, g = 1, a = a, w = w)
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
#' w     <- c(0.6, 1, 1.5, 2, 3.2, 10) ; names(w) <- w
#' olty  <- c(2, 1, 2, 1, 2, 1, 1)
#' olwd  <- c(1, 1, 2, 2, 3, 3, 2)
#' ocol  <- c(2, 2, 4, 4, 3, 3, 1)
#' lleg  <- c("logit(0.999) = 6.9", "logit(0.99)   = 4.6", "logit(0.95)   = 2.9", 
#'            "logit(0.50)   = 0", "logit(0.05)   = -2.9", "logit(0.01)   = -4.6", 
#'            "logit(0.001) = -6.9  ")
#' op    <- par(mfrow=c(2,2), mgp=c(1.5,0.8,0), mar=c(3,3,2,1))
#' 
#' plot(x, plogis(x, scale = 2), type = "n", lwd = 2, ylim = c(0, 1),
#'      xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
#'      main = "pkiener2(q, m, g, a=2, w=...)")
#' for (i in 1:length(w)) lines(x, pkiener2(x, a = 2, w = w[i]), 
#'        lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' legend("topleft", title = expression(omega), legend = c(w), 
#'        cex = 0.7, inset = 0.02, lty = olty, lwd = olwd, col = ocol )
#' 
#' plot(x, dlogis(x, scale = 2), type = "n", lwd = 2, ylim = c(0, 0.17),
#'      xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
#'      main = "dkiener2(q, m, g, a=2, w=...)")
#' for (i in 1:length(w)) lines(x, dkiener2(x, a = 2, w = w[i]), 
#'        lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' legend("topright", title = expression(omega), legend = c(w), 
#'        cex = 0.7, inset = 0.02, lty = olty, lwd = olwd, col = ocol )
#' 
#' plot(x, x/2, type = "n", lwd = 1, ylim = c(-7.5, 7.5), yaxt="n", xaxs = "i", 
#'      yaxs = "i", xlab = "", ylab = "", 
#'      main = "logit(pkiener2(q, m, g, a=2, w=...))")
#' axis(2, las=1, at=c(-6.9, -4.6, -2.9, 0, 2.9, 4.6, 6.9) )
#' for (i in 1:length(w)) lines(x, lkiener2(x, a = 2, w = w[i]),  
#'        lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' legend("topleft", legend = lleg, cex = 0.7, inset = 0.02 )
#' legend("bottomright", title = expression(omega), legend = c(w), 
#'        cex = 0.7, inset = 0.02, lty = c(olty), lwd = c(olwd), col = c(ocol) )
#' 
#' plot(x, dlogis(x, scale = 2, log=TRUE), type = "n", lwd = 2, ylim = c(-8, -1.5),
#'      xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
#'      main = "log(dkiener2(q, m, g, a=2, w=...))")
#' for (i in 1:length(w)) lines(x, dkiener2(x, a = 2, w = w[i], log=TRUE), 
#'        lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' legend("bottom", title = expression(omega), legend = c(w), 
#'        cex = 0.7, inset = 0.02, lty = olty, lwd = olwd, col = ocol )
#' ### End example 3
#' 
#' 
#' ### Example 4 (four plots: quantile, derivate, density and quantiles from p)
#' p     <- ppoints(199, a=0)
#' w     <- c(0.6, 1, 1.5, 2, 3.2, 10) ; names(w) <- w ; w
#' op    <- par(mfrow=c(2,2), mgp=c(1.5,0.8,0), mar=c(3,3,2,1))
#' 
#' plot(p, qlogis(p, scale = 2), type = "l", lwd = 2, xlim = c(0, 1), 
#'      ylim = c(-15, 15), xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
#'      main = "qkiener2(p, m, g, a=2, w=...)")
#' for (i in 1:length(w)) lines(p, qkiener2(p, a = 2, w = w[i]), 
#'           lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' legend("topleft", title = expression(omega), legend = c(w, "qlogis(x/2)"), 
#'           inset = 0.02, lty = olty, lwd = olwd, col = ocol, cex = 0.7 )
#' 
#' plot(p, 2/p/(1-p), type = "l", lwd = 2, xlim = c(0, 1), ylim = c(0, 100),
#'      xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
#'      main = "dqkiener2(p, m, g, a=2, w=...)")
#' for (i in 1:length(w)) lines(p, dqkiener2(p, a = 2, w = w[i]), 
#'           lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' legend("top", title = expression(omega), legend = c(w, "p*(1-p)/2"), 
#'           inset = 0.02, lty = olty, lwd = olwd, col = ocol, cex = 0.7 )
#' 
#' plot(qlogis(p, scale = 2), p*(1-p)/2, type = "l", lwd = 2, xlim = c(-15, 15), 
#'      ylim = c(0, 0.18), xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
#'      main = "qkiener2, dpkiener2(p, m, g, a=2, w=...)")
#' for (i in 1:length(w)) { 
#'      lines(qkiener2(p, a = 2, w = w[i]), dpkiener2(p, a = 2, w = w[i]),
#'            lty = olty[i], lwd = olwd[i], col = ocol[i] ) }
#' legend("topleft", title = expression(omega), legend = c(w, "p*(1-p)/2"), 
#'           inset = 0.02, lty = olty, lwd = olwd, col = ocol, cex = 0.7 )
#' 
#' plot(qlogis(p, scale = 2), p, type = "l", lwd = 2, xlim = c(-15, 15), 
#'      ylim = c(0, 1), xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
#'      main = "inverse axis qkiener2(p, m, g, a=2, w=...)")
#' for (i in 1:length(w)) lines(qkiener2(p, a = 2, w = w[i]), p,
#'           lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' legend("topleft", title = expression(omega), legend = c(w, "qlogis(x/2)"), 
#'           inset = 0.02, lty = olty, lwd = olwd, col = ocol, cex = 0.7 )
#' ### End example 4
#' 
#' 
#' ### Example 5 (q and VaR, ltm, rtm, and ES)
#' pp <- c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 
#'         0.10, 0.20, 0.35, 0.5, 0.65, 0.80, 0.90,
#'         0.95, 0.975, 0.99, 0.995, 0.9975, 0.999)
#' m <- -10 ; g <- 1 ; a <- 5 ; w = 3 
#' k <- aw2k(a, w) ; d <- aw2d(a, w) ; e <- aw2e(a, w)
#' round(c(m = m, g = g, a = a, k = k, w = w, d = d, e = e), 2) 
#' plot(qkiener2(pp, m, g, a, w), pp, type = "b")
#' round(cbind(p = pp, "1-p" = 1-pp,
#' 	q   =   qkiener2(pp, m, g, a, w), 
#' 	ltm = ltmkiener2(pp, m, g, a, w), 
#' 	rtm = rtmkiener2(pp, m, g, a, w), 
#' 	ES  =  eskiener2(pp, m, g, a, w), 
#' 	VaR = varkiener2(pp, m, g, a, w)), 4)
#' round(kmean(c(m, g, a, w), model = "K2"), 4) # limit value for ltm and rtm
#' round(cbind(p = pp, "1-p" = 1-pp, 
#' 	q   =   qkiener2(pp, m, g, a, w, lower.tail = FALSE), 
#' 	ltm = ltmkiener2(pp, m, g, a, w, lower.tail = FALSE), 
#' 	rtm = rtmkiener2(pp, m, g, a, w, lower.tail = FALSE), 
#' 	ES  =  eskiener2(pp, m, g, a, w, lower.tail = FALSE), 
#' 	VaR = varkiener2(pp, m, g, a, w, lower.tail = FALSE)), 4)
#' ### End example 5
#' 
#' 
#' @name kiener2
NULL

#' @export
#' @rdname kiener2
dkiener2 <- function(x, m = 0, g = 1, a = 3.2, w = 3.2, log = FALSE) {
	lp <-  lkiener2(x,  m, g, a, w)
	v  <- dlkiener2(lp, m, g, a, w)
	if(log) return(log(v)) else return(v)
}

#' @export
#' @rdname kiener2
pkiener2 <- function(q, m = 0, g = 1, a = 3.2, w = 3.2, 
                     lower.tail = TRUE, log.p = FALSE) {
	lp <- lkiener2(x = q, m, g, a, w)
	if(lower.tail) v <- invlogit(lp) else v <- 1 - invlogit(lp)
	if(log.p) return(log(v)) else return(v)
}

#' @export
#' @rdname kiener2
qkiener2 <- function(p, m = 0, g = 1, a = 3.2, w = 3.2, 
                     lower.tail = TRUE, log.p = FALSE) {
	if(log.p) p <- exp(p) else p <- p
	if(lower.tail) p <- p else p <- 1-p
	k <- aw2k(a, w)
	v <- m + g * k * (- exp(-logit(p)/a) + exp(logit(p)/w) )
	return(v)
}

#' @export
#' @rdname kiener2
rkiener2 <- function(n, m = 0, g = 1, a = 3.2, w = 3.2) {
	p <- runif(n)
	v <- qkiener2(p, m, g, a, w)
	return(v)
}

#' @export
#' @rdname kiener2
dpkiener2 <- function(p, m = 0, g = 1, a = 3.2, w = 3.2, log = FALSE) {
	k <- aw2k(a, w)
	v <- p * (1 - p) / k / g / ( exp(-logit(p)/a)/a + exp(logit(p)/w)/w )
	if(log) return(log(v)) else return(v)
}

#' @export
#' @rdname kiener2
dqkiener2 <- function(p, m = 0, g = 1, a = 3.2, w = 3.2, log = FALSE) {
# Compute dX/dp
	k <- aw2k(a, w)
	v <- k * g / p / (1 - p) * ( exp(-logit(p)/a)/a + exp(logit(p)/w)/w )
	if(log) return(log(v)) else return(v)
}


#' @export
#' @rdname kiener2
lkiener2 <- function(x, m = 0, g = 1, a = 3.2, w = 3.2) { 
	k      <- min(a, w)
	lp.ini <- lkiener1(x, m, g, k)
	f      <- function(lp) sum( ( x - qlkiener2(lp, m, g, a, w) )^2 )
	lp.fin <- nlm(f, lp.ini)
	v      <- lp.fin$estimate
	return(v)
}

#' @export
#' @rdname kiener2
dlkiener2 <- function(lp, m = 0, g = 1, a = 3.2, w = 3.2, log = FALSE) {
	p <- invlogit(lp)
	k <- aw2k(a, w)
	v <- p * (1 - p) / k / g / ( exp(-lp/a)/a + exp(lp/w)/w )
	if(log) return(log(v)) else return(v)
}

#' @export
#' @rdname kiener2
qlkiener2 <- function(lp, m = 0, g = 1, a = 3.2, w = 3.2, lower.tail = TRUE ) {
	if(lower.tail) lp <- lp else lp <- -lp
	k <- aw2k(a, w)
	v <- m + g * k * ( - exp(-lp/a) + exp(lp/w) )
	return(v)
}

#' @export
#' @rdname kiener2
varkiener2 <- function(p, m = 0, g = 1, a = 3.2, w = 3.2, 
                      lower.tail = TRUE, log.p = FALSE) {
	p   <- if(log.p) {exp(p)} else {p}
	p   <- if(lower.tail) {p} else {1-p}
	va  <- p
	for (i in seq_along(p)) {
		va[i] <- ifelse(p[i] <= 0.5, 
					- qkiener2(p[i], m, g, a, w),
					  qkiener2(p[i], m, g, a, w))
	}	
return(va)
}

#' @export
#' @rdname kiener2
ltmkiener2 <- function(p, m = 0, g = 1, a = 3.2, w = 3.2, 
                       lower.tail = TRUE, log.p = FALSE) {
	p  <- if (log.p) {exp(p)} else {p}
	k  <- aw2k(a, w)
	ltm <- if (lower.tail) {
		m+g*k/p*(
			-pbeta(p, 1-1/a, 1+1/a)*beta(1-1/a, 1+1/a)
			+pbeta(p, 1+1/w, 1-1/w)*beta(1+1/w, 1-1/w))	
	} else {
		m+g*k/p*(
			-pbeta(p, 1+1/a, 1-1/a)*beta(1+1/a, 1-1/a)
			+pbeta(p, 1-1/w, 1+1/w)*beta(1-1/w, 1+1/w))
	}
return(ltm)
}

#' @export
#' @rdname kiener2
rtmkiener2 <- function(p, m = 0, g = 1, a = 3.2, w = 3.2, 
                       lower.tail = TRUE, log.p = FALSE) {
	p  <- if (log.p) {exp(p)} else {p}
	k  <- aw2k(a, w)
	rtm <- if (!lower.tail) {
		m+g*k/(1-p)*(
			-pbeta(1-p, 1-1/a, 1+1/a)*beta(1-1/a, 1+1/a)
			+pbeta(1-p, 1+1/w, 1-1/w)*beta(1+1/w, 1-1/w))	
	} else {
		m+g*k/(1-p)*(
			-pbeta(1-p, 1+1/a, 1-1/a)*beta(1+1/a, 1-1/a)
			+pbeta(1-p, 1-1/w, 1+1/w)*beta(1-1/w, 1+1/w))
	}
return(rtm)
}

#' @export
#' @rdname kiener2
dtmqkiener2 <- function(p, m = 0, g = 1, a = 3.2, w = 3.2, 
                      lower.tail = TRUE, log.p = FALSE) {
	dtmq <- p
	for (i in seq_along(p)) {
		dtmq[i] <- ifelse(p[i] <= 0.5, 
			ltmkiener2(p[i], m, g, a, w, lower.tail, log.p) 
			- qkiener2(p[i], m, g, a, w, lower.tail, log.p),
			rtmkiener2(p[i], m, g, a, w, lower.tail, log.p) 
			- qkiener2(p[i], m, g, a, w, lower.tail, log.p))	
	}
return(dtmq)
}

#' @export
#' @rdname kiener2
eskiener2 <- function(p, m = 0, g = 1, a = 3.2, w = 3.2, 
                      lower.tail = TRUE, log.p = FALSE) {
	p   <- if (log.p) {exp(p)} else {p}
	p   <- if(lower.tail) {p} else {1-p}
	es  <- p
	for (i in seq_along(p)) {
		es[i] <- ifelse(p[i] <= 0.5, 
					- ltmkiener2(p[i], m, g, a, w),
					  rtmkiener2(p[i], m, g, a, w))
	}
return(es)
}



