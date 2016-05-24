

#' @include e_conversion.R



#' @title Symmetric Kiener Distribution K1
#'
#' @description
#' Density, distribution function, quantile function, random generation
#' and additional formulae for symmetric Kiener distribution K1. 
#' This distribution is similar to the power hyperbola logistic distribution 
#' but with additional parameters for location (\code{m}) and scale (\code{g}).
#'
#' @param    x    vector of quantiles.
#' @param    q	  vector of quantiles.
#' @param    m	  numeric. The median.
#' @param    g	  numeric. The scale parameter, preferably strictly positive.
#' @param    k	  numeric. The tail parameter, preferably strictly positive.
#' @param    p	  vector of probabilities.
#' @param    lp	  vector of logit of probabilities.
#' @param    n	  number of observations. If length(n) > 1, the length is  
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
#' Kiener distributions \code{K1(m, g, k, ...)} describe distributions  
#' with symmetric left and right fat tails with tail parameter \code{k}. 
#' This parameter is the power exponent mentionned in Pareto formula and 
#' Karamata theorems.
#' 
#' \code{m} is the median of the distribution. \code{g} is the scale parameter 
#' and the inverse of the density at the median: \eqn{ g = 1 / 8 / f(m) }.
#' As a first estimate, it is approximatively one fourth of the standard 
#' deviation \eqn{ g  \approx \sigma / 4 } but is independant from it.
#'
#' \code{dkiener1} function is defined for x in (-Inf, +Inf) by: 
#'      \deqn{ dkiener1(x, m, g, k) = 
#'                 1 / 4 / g / cosh(  ashp((x - m)/g, k) ) 
#'                      / (1 + cosh( kashp((x - m)/g, k))) }
#'
#' \code{pkiener1} function is defined for q in (-Inf, +Inf) by: 
#'       \deqn{ pkiener1(q, m, g, k) = 1/(1 + exp(- kashp((q - m)/g, k))) }
#'
#' \code{qkiener1} function is defined for p in (0, 1) by: 
#'       \deqn{ qkiener1(p, m, g, k) = m + 2 * g * k * sinh( logit(p)/k ) }
#'
#' \code{rkiener1} generates \code{n} random quantiles.
#'
#' In addition to the classical d, p, q, r functions, the prefixes 
#' dp, dq, l, dl, ql are also provided.
#'
#' \code{dpkiener1} is the density function calculated from the probability p. 
#' It is defined for p in (0, 1) by: 
#'    \deqn{ dpkiener1(p, m, g, k) = p * (1 - p) / 2 / g / cosh( logit(p)/k ) }
#'
#' \code{dqkiener1} is the derivate of the quantile function calculated from 
#' the probability p. It is defined for p in (0, 1) by: 
#'    \deqn{ dqkiener1(p, m, g, k) = 2 * g / p / (1 - p) * cosh( logit(p)/k ) }
#'
#' \code{lkiener1} function is equivalent to kashp function but with additional 
#' parameters \code{m} and \code{g}. Being computed from the x (or q) vector, 
#' it can be compared to the logit of the empirical probability logit(p) 
#' through a nonlinear regression with ordinary or weighted least squares 
#' to estimate the distribution parameters. 
#' It is defined for x in (-Inf, +Inf) by:
#'    \deqn{ lkiener1(x, m, g, k) = kashp((x - m)/g, k) }
#'
#' \code{dlkiener1} is the density function calculated from the logit of the 
#' probability lp = logit(p). It is defined for lp in (-Inf, +Inf) by: 
#'    \deqn{ dlkiener1(lp, m, g, k) = p * (1 - p) / 2 / g / cosh( lp/k ) }
#'
#' \code{qlkiener1} is the quantile function calculated from the logit of the 
#' probability lp = logit(p). It is defined for lp in (-Inf, +Inf) by: 
#'    \deqn{ qlkiener1(lp, m, g, k) = m + g * k * 2 * sinh( lp/k ) }
#' 
#' \code{varkiener1} designates the Value a-risk and turns negative numbers 
#' into positive numbers with the following rule:
#'    \deqn{ varkiener1 <- if(p <= 0.5) (- qkiener1) else (qkiener1) }
#' Usual values in finance are \code{p = 0.01}, \code{p = 0.05}, \code{p = 0.95} and 
#' \code{p = 0.99}. \code{lower.tail = FALSE} uses \code{1-p} rather than {p}.
#'
#' \code{ltmkiener1}, \code{rtmkiener1} and \code{eskiener1} are respectively the 
#' left tail mean, the right tail mean and the expected shortfall of the distribution 
#' (sometimes called average VaR, conditional VaR or tail VaR). 
#' Left tail mean is the integrale from \code{-Inf} to \code{p} of the quantile function 
#' \code{qkiener1} divided by \code{p}.
#' Right tail mean is the integrale from \code{p} to \code{+Inf} of the quantile function 
#' \code{qkiener1} divided by 1-p.
#' Expected shortfall turns negative numbers into positive numbers with the following rule:
#'    \deqn{ eskiener1 <- if(p <= 0.5) (- ltmkiener1) else (rtmkiener1) }
#' Usual values in finance are \code{p = 0.01}, \code{p = 0.025}, \code{p = 0.975} and 
#' \code{p = 0.99}. \code{lower.tail = FALSE} uses \code{1-p} rather than {p}.
#'
#' \code{dtmqkiener1} is the difference between the left tail mean and the quantile 
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
#' Power hyperbola logistic distribution \code{\link{logishp}}, 
#' asymmetric Kiener distributions K2, K3 and K4  
#' \code{\link{kiener2}}, \code{\link{kiener3}}, \code{\link{kiener4}}, 
#' regression function \code{\link{regkienerLX}}.
#'
#' @examples
#' 
#' require(graphics)
#' 
#' ### Example 1
#' pp <- c(ppoints(11, a = 1), NA, NaN) ; pp
#' qkiener1(p = pp, k = 4)
#' 
#' 
#' ### Example 2: Try various value of k = 1.5, 3, 5, 10
#' k       <- 5  # 1.5, 3, 5, 10
#' set.seed(2014)
#' mainTC  <- paste("qkiener1(p, m = 0, g = 1, k = ", k, ")")
#' mainsum <- paste("cumulated qkiener1(p, m = 0, g = 1, k = ", k, ")")
#' T       <- 500
#' C       <- 4
#' TC      <- qkiener1(p = runif(T*C), m = 0, g = 1, k = k)
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
#' x  <- q  <- seq(-15, 15, length.out=101)
#' k     <- c(0.6, 1, 1.5, 2, 3.2, 10) ; names(k) <- k ; k
#' olty  <- c(2, 1, 2, 1, 2, 1, 1)
#' olwd  <- c(1, 1, 2, 2, 3, 3, 2)
#' ocol  <- c(2, 2, 4, 4, 3, 3, 1)
#' lleg  <- c("logit(0.999) = 6.9", "logit(0.99)   = 4.6", "logit(0.95)   = 2.9", 
#'            "logit(0.50)   = 0", "logit(0.05)   = -2.9", "logit(0.01)   = -4.6", 
#'            "logit(0.001) = -6.9  ")
#' op    <- par(mfrow=c(2,2), mgp=c(1.5,0.8,0), mar=c(3,3,2,1))
#' 
#' plot(x, plogis(x, scale = 2), type = "b", lwd = 2, ylim = c(0, 1),
#'      xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
#'      main = "pkiener1(q, m, g, k)")
#' for (i in 1:length(k)) lines(x, pkiener1(x, k = k[i]), 
#'        lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' legend("topleft", title = expression(kappa), legend = c(k, "logistic"), 
#'        cex = 0.7, inset = 0.02, lty = olty, lwd = olwd, col = ocol )
#' 
#' plot(x, dlogis(x, scale = 2), type = "b", lwd = 2, ylim = c(0, 0.14), 
#'      xaxs = "i", yaxs = "i", xlab = "", ylab = "", main = "dkiener1(x, m, g, k)" )
#' for (i in 1:length(k)) lines(x, dkiener1(x, k = k[i]), 
#'        lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' legend("topright", title = expression(kappa), legend = c(k, "logistic"), 
#'        cex = 0.7, inset = 0.02, lty = olty, lwd = olwd, col = ocol )
#' 
#' plot(x, x/2, type = "b", lwd = 2, ylim = c(-7.5, 7.5), yaxt="n", xaxs = "i", 
#'      yaxs = "i", xlab = "", ylab = "", main = "logit(pkiener1(q, m, g, k))")
#' axis(2, las=1, at=c(-6.9, -4.6, -2.9, 0, 2.9, 4.6, 6.9) )
#' for (i in 1:length(k)) lines(x, lkiener1(x, k = k[i]),  
#'        lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' lines(x, logit(pnorm(x, 0, 3.192)), type="l", lty=1, lwd=3, col=7) # erfx
#' legend("topleft", legend = lleg, cex = 0.7, inset = 0.02 )
#' legend("bottomright", title = expression(kappa), 
#'        legend = c(k, "logistic", "Gauss"), cex = 0.7, inset = 0.02, 
#'        lty = c(olty, 1), lwd = c(olwd, 3), col = c(ocol , 7) )
#' 
#' plot(x, log(dlogis(x, scale = 2)), lwd = 2, type = "b", ylim = c(-8, -1.5), 
#'      xaxs = "i", yaxs = "i", xlab = "", ylab = "", main = "log(dkiener1(x, m, g, k))") 
#' for (i in 1:length(k)) lines(x, log(dkiener1(x, k = k[i])),  
#'        lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' lines(x, dnorm(x, 0, 3.192, log = TRUE), type = "l", lty = 1, lwd = 3, col = 7)
#' legend("bottom", title = expression(kappa), legend = c(k, "logistic", "Gauss"), 
#'        cex = 0.7, inset = 0.02, lty = c(olty, 1), lwd = c(olwd, 3), col = c(ocol , 7) )
#' ### End example 3
#' 
#' 
#' ### Example 4 (four plots: quantile, derivate, density and quantiles from p)
#' p   <- ppoints(199, a=0)
#' k   <- c(0.6, 1, 1.5, 2, 3.2, 10) ; names(k) <- k ; k
#' op  <- par(mfrow=c(2,2), mgp=c(1.5,0.8,0), mar=c(3,3,2,1))
#' plot(p, qlogis(p, scale = 2), type = "o", lwd = 2, ylim = c(-15, 15),
#'      xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
#'      main = "qkiener1(p, m, g, k)")
#' for (i in 1:length(k)) lines(p, qkiener1(p, k = k[i]), 
#'           lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' legend("topleft", title = expression(kappa), legend = c(k, "qlogis(x/2)"), 
#'           inset = 0.02, lty = olty, lwd = olwd, col = ocol, cex = 0.7 )
#' 
#' plot(p, 2/p/(1-p), type = "o", lwd = 2, xlim = c(0, 1), ylim = c(0, 100),
#'      xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
#'      main = "dqkiener1(p, m, g, k)")
#' for (i in 1:length(k)) lines(p, dqkiener1(p, k = k[i]), 
#'           lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' legend("top", title = expression(kappa), legend = c(k, "p*(1-p)/2"), 
#'           inset = 0.02, lty = olty, lwd = olwd, col = ocol, cex = 0.7 )
#' 
#' plot(qlogis(p, scale = 2), p*(1-p)/2, type = "o", lwd = 2, xlim = c(-15, 15), 
#'      ylim = c(0, 0.14), xaxs = "i", yaxs = "i", xlab = "", ylab = "", 
#'      main = "qkiener1, dpkiener1(p, m, g, k)")
#' for (i in 1:length(k)) lines(qkiener1(p, k = k[i]), dpkiener1(p, k = k[i]),
#'           lty = olty[i], lwd = olwd[i], col = ocol[i] )
#' legend("topleft", title = expression(kappa), legend = c(k, "p*(1-p)/2"), 
#'           inset = 0.02, lty = olty, lwd = olwd, col = ocol, cex = 0.7 )
#' ### End example 4
#' 
#' 
#' ### Example 5 (q and VaR, ltm, rtm, and ES)
#' pp <- c(0.001, 0.0025, 0.005, 0.01, 0.025, 0.05, 
#'         0.10, 0.20, 0.35, 0.5, 0.65, 0.80, 0.90,
#'         0.95, 0.975, 0.99, 0.995, 0.9975, 0.999)
#' m <- -10 ; g <- 1 ; k <- 4
#' round(c(m = m, g = g, a = k, k = k, w = k, d = 0, e = 0), 2) 
#' plot(qkiener1(pp, m, g, k), pp, type = "b")
#' round(cbind(p = pp, "1-p" = 1-pp, 
#' 	q   =   qkiener1(pp, m, g, k), 
#' 	ltm = ltmkiener1(pp, m, g, k), 
#' 	rtm = rtmkiener1(pp, m, g, k), 
#' 	es  =  eskiener1(pp, m, g, k), 
#' 	VaR = varkiener1(pp, m, g, k)), 4)
#' round(kmean(c(m, g, k), model = "K1"), 4) # limit value of ltm, rtm
#' round(cbind(p = pp, "1-p" = 1-pp, 
#' 	q   =   qkiener1(pp, m, g, k, lower.tail = FALSE), 
#' 	ltm = ltmkiener1(pp, m, g, k, lower.tail = FALSE), 
#' 	rtm = rtmkiener1(pp, m, g, k, lower.tail = FALSE), 
#' 	es  =  eskiener1(pp, m, g, k, lower.tail = FALSE), 
#' 	VaR = varkiener1(pp, m, g, k, lower.tail = FALSE)), 4)
#' ### End example 5
#' 
#' 
#' @name kiener1
NULL

#' @export
#' @rdname kiener1
dkiener1 <- function(x, m = 0, g = 1, k = 3.2, log = FALSE) {
	ash <- asinh((x - m) / g / 2 / k)
	v <- 1 / 4 / g / cosh(ash) / (1 + cosh(k * ash))
	if(log) return(log(v)) else return(v)
}

#' @export
#' @rdname kiener1
pkiener1 <- function(q, m = 0, g = 1, k = 3.2, lower.tail = TRUE, log.p = FALSE) {
	v <-  1/(1 + exp(- kashp((q - m)/g, k)))
	if(lower.tail) v <- v else v <- 1-v
	if(log.p) return(log(v)) else return(v)
}

#' @export
#' @rdname kiener1
qkiener1 <- function(p, m = 0, g = 1, k = 3.2, lower.tail = TRUE, log.p = FALSE) {
	if(log.p) p <- exp(p) else p <- p
	if(lower.tail) p <- p else p <- 1-p
	v <- m + g * k * 2 * sinh( logit(p)/k )
	return(v)
}

#' @export
#' @rdname kiener1
rkiener1 <- function(n, m = 0, g = 1, k = 3.2) {
	p <- runif(n)
	v <- qkiener1(p, m, g, k, lower.tail = TRUE, log.p = FALSE)
	return(v)
}

#' @export
#' @rdname kiener1
dpkiener1 <- function(p, m = 0, g = 1, k = 3.2, log = FALSE) { 
	v <- p * (1 - p) / 2 / g / cosh( logit(p)/k )
	if(log) return(log(v)) else return(v)
}

#' @export
#' @rdname kiener1
dqkiener1 <- function(p, m = 0, g = 1, k = 3.2, log = FALSE) { 
	v <- 2 * g / p / (1 - p) * cosh( logit(p)/k )
	if(log) return(log(v)) else return(v)
}

#' @export
#' @rdname kiener1
lkiener1 <- function(x, m = 0, g = 1, k = 3.2) { 
	kashp((x - m)/g, k) 
}

#' @export
#' @rdname kiener1
dlkiener1 <- function(lp, m = 0, g = 1, k = 3.2, log = FALSE) {
	p <- plogis(lp)       # = invlogit
	v <- p * (1 - p) / 2 / g / cosh( lp/k )
	if(log) return(log(v)) else return(v)
}

#' @export
#' @rdname kiener1
qlkiener1 <- function(lp, m = 0, g = 1, k = 3.2, lower.tail = TRUE ) {
	if(lower.tail) lp <- lp else lp <- -lp
	v <- m + g * k * 2 * sinh( lp/k )
	return(v)
}

#' @export
#' @rdname kiener1
varkiener1 <- function(p, m = 0, g = 1, k = 3.2, 
                      lower.tail = TRUE, log.p = FALSE) {
	p   <- if(log.p) {exp(p)} else {p}
	p   <- if(lower.tail) {p} else {1-p}
	va  <- p
	for (i in seq_along(p)) {
		va[i] <- ifelse(p[i] <= 0.5, 
					- qkiener1(p[i], m, g, k),
					  qkiener1(p[i], m, g, k))
	}
return(va)
}

#' @export
#' @rdname kiener1
ltmkiener1 <- function(p, m = 0, g = 1, k = 3.2, 
                       lower.tail = TRUE, log.p = FALSE) {
	p   <- if (log.p) {exp(p)} else {p}
	ltm <- if (lower.tail) {
		m+g*k/p*(
			-pbeta(p, 1-1/k, 1+1/k)*beta(1-1/k, 1+1/k)
			+pbeta(p, 1+1/k, 1-1/k)*beta(1+1/k, 1-1/k))	
	} else {
		m+g*k/p*(
			-pbeta(p, 1+1/k, 1-1/k)*beta(1+1/k, 1-1/k)
			+pbeta(p, 1-1/k, 1+1/k)*beta(1-1/k, 1+1/k))
	}
return(ltm)
}

#' @export
#' @rdname kiener1
rtmkiener1 <- function(p, m = 0, g = 1, k = 3.2, 
                       lower.tail = TRUE, log.p = FALSE) {
	p   <- if(log.p) {exp(p)} else {p}
	rtm <- if (!lower.tail) {
		m+g*k/(1-p)*(
			-pbeta(1-p, 1-1/k, 1+1/k)*beta(1-1/k, 1+1/k)
			+pbeta(1-p, 1+1/k, 1-1/k)*beta(1+1/k, 1-1/k))	
	} else {
		m+g*k/(1-p)*(
			-pbeta(1-p, 1+1/k, 1-1/k)*beta(1+1/k, 1-1/k)
			+pbeta(1-p, 1-1/k, 1+1/k)*beta(1-1/k, 1+1/k))
	}
return(rtm)
}

#' @export
#' @rdname kiener1
dtmqkiener1 <- function(p, m = 0, g = 1, k = 3.2, 
                        lower.tail = TRUE, log.p = FALSE) {
	dtmq <- p
	for (i in seq_along(p)) {
		dtmq[i] <- ifelse(p[i] <= 0.5, 
			ltmkiener1(p[i], m, g, k, lower.tail, log.p) 
			- qkiener1(p[i], m, g, k, lower.tail, log.p),
			rtmkiener1(p[i], m, g, k, lower.tail, log.p) 
			- qkiener1(p[i], m, g, k, lower.tail, log.p))	
	}
return(dtmq)
}

#' @export
#' @rdname kiener1
eskiener1 <- function(p, m = 0, g = 1, k = 3.2, 
                      lower.tail = TRUE, log.p = FALSE) {
	p   <- if (log.p) {exp(p)} else {p}
	p   <- if(lower.tail) {p} else {1-p}
	es  <- p
	for (i in seq_along(p)) {
		es[i] <- ifelse(p[i] <= 0.5, 
					- ltmkiener1(p[i], m, g, k),
					  rtmkiener1(p[i], m, g, k))
	}
return(es)
}





