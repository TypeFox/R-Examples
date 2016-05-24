# Copyright 2015-2015 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav
#
# This file is part of PDQutils.
#
# PDQutils is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PDQutils is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with PDQutils.  If not, see <http://www.gnu.org/licenses/>.

# Created: 2015.02.07
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav
# Comments: Steven E. Pav

# for the Hermite Polynomials
require(orthopolynom)
require(moments)


# suppose raw.moments[k] is the kth raw moment of X;
# here we will compute the kth raw moment of X+del.
# n.b.
# E[(x+del)^k] = E[x^k + del choose(k,1) x^{k-1} + ... del^k ]
.shift_moments <- function(raw.moments,del) {
	nmom <- length(raw.moments)-1
	shf.moments <- raw.moments
	shf.moments[2] <- raw.moments[2] + del
	for (k in 2:nmom) {
		tot <- 0
		for (j in 0:k) {
			tot <- tot + choose(k, j) * del^(k-j) * raw.moments[j+1]
		}
		shf.moments[k+1] <- tot
	}
	return(shf.moments)
}
# suppose raw.moments[k] is the kth raw moment of X;
# here we will compute the kth raw moment of a * X.
# n.b.
# E[(ax)^k] = a^k E[x^k]
.scale_moments <- function(raw.moments,k) {
	nmom <- length(raw.moments)-1
	scl.moments <- raw.moments * (k^(0:nmom))
	return(scl.moments)
}

.gca_setup <- function(x,raw.moments,support=NULL, basis=c('normal','gamma','beta','arcsine','wigner'), basepar=NULL) {
	basis <- tolower(match.arg(basis))
	# the zeroth moment
	raw.moments <- c(1,raw.moments)

	# guess support:#FOLDUP
	if (is.null(support)) { 
		support <- switch(basis,
											"normal"=c(-Inf,Inf),
											"gamma"=c(0,Inf),
											"beta"=c(0,1),
											"arcsine"=c(-1,1),
											"wigner"=c(-1,1))
	}#UNFOLD
	support <- sort(support)
	# make these special cases of beta:#FOLDUP
	if (basis == 'arcsine') {
		basepar = list(shape1=0.5,shape2=0.5)
		basis = 'beta'
	} else if (basis == 'wigner') {
		basepar = list(shape1=1.5,shape2=1.5)
		basis = 'beta'
	}#UNFOLD
	# shift, scale X, modify the moments, compute final scaling factor#FOLDUP
	if (basis == 'normal') {
		mu <- raw.moments[2]
		sigma <- sqrt(raw.moments[3] - mu^2)
		x <- (x - mu) / sigma
		moments <- .shift_moments(raw.moments,-mu)
		moments <- .scale_moments(moments,1/sigma)
		scalby <- 1/sigma
		support <- (support - mu)/sigma
	} else if (basis == 'gamma') {
		llim = min(support)
		x <- (x - llim)
		moments <- .shift_moments(raw.moments,-llim)
		support <- (support - llim)
		scalby <- 1
	} else if (basis == 'beta') {
		ulim = max(support)
		llim = min(support)
		mu <- 0.5 * (ulim + llim)
		sigma <- 0.5 * (ulim - llim)
		x <- (x - mu) / sigma
		moments <- .shift_moments(raw.moments,-mu)
		moments <- .scale_moments(moments,1/sigma)
		scalby <- 1/sigma
		support <- c(-1,1)
	} else { stop('badCode') }#UNFOLD  # nocov
	# guess the base distribution parameters, from the moments?#FOLDUP
	if (is.null(basepar)) {
		if (basis == 'gamma') {
			# first two uncentered moments for gamma are k theta and k theta^2 + (k theta)^2
			theta <- (moments[3]/moments[2]) - moments[2]
			k <- moments[2] / theta
			basepar <- list(shape=k,scale=theta)
		} else if (basis == 'beta') {
			# compute a and b
			mu <- moments[2]
			s2 <- moments[3] - moments[2]^2
			# shift back to [0,1]
			mu <- (mu + 1) / 2 
			s2 <- s2 / 4
			# second moment
			mu2 <- s2 + mu^2
			# solve for b, a
			b <- (mu - mu2) * (1 - mu) / s2
			a <- b * mu / (1-mu)
			# n.b. the reverse
			basepar <- list(shape2=b,shape1=a)
		}
	}#UNFOLD
	# rescale gammas#FOLDUP
	if (basis == 'gamma') {
		x <- x / basepar$scale
		moments <- .scale_moments(moments,1/basepar$scale)
		scalby <- scalby / basepar$scale
		support <- support / basepar$scale
		basepar$scale <- 1
	}#UNFOLD

	order.max <- length(moments)-1
	orders <- seq(0,order.max)
	
	if (basis == 'normal') {
		wt <- dnorm
		# the orthogonal polynomials
		poly <- orthopolynom::hermite.he.polynomials(order.max, normalized=FALSE)
		hn <- factorial(orders)
		intpoly <- c(function(y) { as.numeric(poly[[1]]) * pnorm(y) },
								 lapply(poly[1:(order.max)],function(pol) { function(y) { -dnorm(y) * as.function(pol)(y) } }) )
	} else if (basis == 'gamma') {
		alpha <- basepar$shape - 1
		wt <- function(x) { dgamma(x,shape=alpha+1,scale=1) }
		poly <- orthopolynom::glaguerre.polynomials(order.max, alpha, normalized=FALSE)
		hn <- exp(lgamma(alpha + 1 + orders) - lgamma(alpha+1) - lfactorial(orders))
		ipoly <- orthopolynom::glaguerre.polynomials(order.max-1, alpha+1, normalized=FALSE)
		intpoly <- c(function(y) { as.numeric(poly[[1]]) * pgamma(y,shape=alpha+1,scale=1) },
								 lapply(1:(order.max),
												function(idx) { function(y) { ((alpha+1)/idx) * dgamma(y,shape=alpha+2,scale=1) * as.function(ipoly[[idx]])(y) } }) )

	} else if (basis == 'beta') {
		palpha <- basepar$shape2 - 1
		pbeta <- basepar$shape1 - 1
		wt <- function(x) { 0.5 * dbeta(0.5 * (x+1),shape2=palpha+1,shape1=pbeta+1) }
	
		poly <- orthopolynom::jacobi.p.polynomials(order.max, alpha=palpha, beta=pbeta, normalized=FALSE)
		hn <- exp(lgamma(orders + palpha + 1) + lgamma(orders + pbeta + 1) - lfactorial(orders) - lgamma(orders + palpha + pbeta + 1) - 
							lbeta(palpha+1,pbeta+1) - log(2*orders+palpha+pbeta+1))

		ipoly <- orthopolynom::jacobi.p.polynomials(order.max-1, alpha=palpha+1, beta=pbeta+1, normalized=FALSE)
		intpoly <- c(function(y) { as.numeric(poly[[1]]) * pbeta(0.5 * (y+1),shape2=palpha+1,shape1=pbeta+1) },
								 lapply(1:(order.max),
												function(idx) { function(y) { 
													(-2/idx) * exp(lbeta(palpha+2,pbeta+2) - lbeta(palpha+1,pbeta+1)) *
													(0.5 * dbeta(0.5 * (x+1),shape1=palpha+2,shape2=pbeta+2)) *
													as.function(ipoly[[idx]])(y) } }))
	} else { stop(paste('badCode: distribution',basis,'unknown')) } # nocov

	retval <- list(x=x,full_moments=moments,support=support,scalby=scalby,
								 order.max=order.max,orders=orders,
								 wt=wt,poly=poly,hn=hn,intpoly=intpoly)
}


#' @title Approximate density and distribution via Gram-Charlier A expansion.
#'
#' @description 
#'
#' Approximate the probability density or cumulative distribution function of a distribution via its raw moments.
#'
#' @template details-gca
#'
#' @usage
#'
#' dapx_gca(x, raw.moments, support=NULL, 
#'  basis=c('normal','gamma','beta','arcsine','wigner'), 
#'  basepar=NULL, log=FALSE)
#'
#' papx_gca(q, raw.moments, support=NULL, 
#'  basis=c('normal','gamma','beta','arcsine','wigner'), 
#'  basepar=NULL, lower.tail=TRUE, log.p=FALSE)
#'
#' @param x where to evaluate the approximate density.
#' @param q where to evaluate the approximate distribution.
#' @param raw.moments an atomic array of the 1st through kth raw moments
#' of the probability distribution. 
#' @param support the support of the density function. It is assumed
#' that the density is zero on the complement of this open interval.
#' This defaults to \code{c(-Inf,Inf)} for the normal basis,
#' \code{c(0,Inf)} for the gamma basis, and
#' \code{c(0,1)} for the Beta, and 
#' \code{c(-1,1)} for the arcsine and wigner.
#' @param basis the basis under which to perform the approximation. \code{'normal'}
#' gives the classical 'A' series expansion around the PDF and CDF of the normal
#' distribution via Hermite polynomials. \code{'gamma'} expands around a
#' gamma distribution with parameters \code{basepar$shape} and
#' \code{basepar$scale}. 
#' \code{'beta'} expands around a beta distribution with parameters
#' \code{basepar$shape1} and \code{basepar$shape2}. 
#' @param basepar the parameters for the base distribution approximation. 
#' If \code{NULL}, the shape and rate are inferred from the first two moments
#' and/or from the \code{support} as appropriate.
#' @param log logical; if TRUE, densities \eqn{f} are given 
#'  as \eqn{\mbox{log}(f)}{log(f)}.
#' @param log.p logical; if TRUE, probabilities p are given 
#'  as \eqn{\mbox{log}(p)}{log(p)}.
#' @param lower.tail whether to compute the lower tail. If false, we approximate the survival function.
#' @return The approximate density at \code{x}, or the approximate CDF at
#' @keywords distribution 
#' @seealso \code{\link{qapx_cf}}
#' @export 
#' @template ref-Jaschke
#' @template ref-Blinnikov
#' @aliases papx_gca 
#' @note 
#'
#' Monotonicity of the CDF is not guaranteed.
#'
#' @examples 
#' # normal distribution:
#' xvals <- seq(-2,2,length.out=501)
#' d1 <- dapx_gca(xvals, c(0,1,0,3,0), basis='normal')
#' d2 <- dnorm(xvals)
#' # they should match:
#' d1 - d2
#'
#' qvals <- seq(-2,2,length.out=501)
#' p1 <- papx_gca(qvals, c(0,1,0,3,0))
#' p2 <- pnorm(qvals)
#' p1 - p2
#'
#' xvals <- seq(-6,6,length.out=501)
#' mu <- 2
#' sigma <- 3
#' raw.moments <- c(2,13,62,475,3182)
#' d1 <- dapx_gca(xvals, raw.moments, basis='normal')
#' d2 <- dnorm(xvals,mean=mu,sd=sigma)
#' \dontrun{
#' plot(xvals,d1)
#' lines(xvals,d2,col='red')
#' }
#' p1 <- papx_gca(xvals, raw.moments, basis='normal')
#' p2 <- pnorm(xvals,mean=mu,sd=sigma)
#' \dontrun{
#' plot(xvals,p1)
#' lines(xvals,p2,col='red')
#' }
#'
#' # for a one-sided distribution, like the chi-square
#' chidf <- 30
#' ords <- seq(1,9)
#' raw.moments <- exp(ords * log(2) + lgamma((chidf/2) + ords) - lgamma(chidf/2))
#' xvals <- seq(0.3,10,length.out=501)
#' d1g <- dapx_gca(xvals, raw.moments, support=c(0,Inf), basis='gamma')
#' d2 <- dchisq(xvals,df=chidf)
#' \dontrun{
#' plot(xvals,d1g)
#' lines(xvals,d2,col='red')
#' }
#' 
#' p1g <- papx_gca(xvals, raw.moments, support=c(0,Inf), basis='gamma')
#' p2 <- pchisq(xvals,df=chidf)
#' \dontrun{
#' plot(xvals,p1g)
#' lines(xvals,p2,col='red')
#' }
#'
#' # for a one-sided distribution, like the log-normal
#' mu <- 2
#' sigma <- 1
#' ords <- seq(1,8)
#' raw.moments <- exp(ords * mu + 0.5 * (sigma*ords)^2)
#' xvals <- seq(0.5,10,length.out=501)
#' d1g <- dapx_gca(xvals, raw.moments, support=c(0,Inf), basis='gamma')
#' d2 <- dnorm(log(xvals),mean=mu,sd=sigma) / xvals
#' \dontrun{
#' 	plot(xvals,d1g)
#' 	lines(xvals,d2,col='red')
#' }
#' @template etc
dapx_gca <- function(x,raw.moments,support=NULL,basis=c('normal','gamma','beta','arcsine','wigner'),basepar=NULL,
										 log=FALSE) {#FOLDUP

	basis <- tolower(match.arg(basis))
	gca <- .gca_setup(x,raw.moments,support,basis,basepar)

	wx <- gca$wt(gca$x)
	retval <- (as.numeric(gca$poly[[1]]) / gca$hn[1]) * wx
	for (iii in c(1:gca$order.max)) {
		ci <- (sum(coef(gca$poly[[iii+1]]) * gca$full_moments[1:(iii+1)])) / gca$hn[iii+1] 
		retval <- retval + ci * wx * (as.function(gca$poly[[iii+1]])(gca$x))
	}
	# adjust back from standardized
	retval <- retval * gca$scalby

	# sanity check; shall I throw a warning?
	retval <- pmax(0,retval)

	# support support
	if (is.finite(min(gca$support))) {
		retval[gca$x <= min(gca$support)] <- 0
	}
	if (is.finite(max(gca$support))) {
		retval[gca$x >= max(gca$support)] <- 0
	}

	# must be a better way to do this ... 
	if (log)
		retval <- log(retval)
	return(retval)
}#UNFOLD
#' @export 
papx_gca <- function(q,raw.moments,support=NULL,basis=c('normal','gamma','beta','arcsine','wigner'),basepar=NULL,
										 lower.tail=TRUE,log.p=FALSE) {#FOLDUP
	basis <- tolower(match.arg(basis))
	gca <- .gca_setup(q,raw.moments,support,basis,basepar)

	retval <- 0
	for (iii in c(0:gca$order.max)) {
		ci <- (sum(coef(gca$poly[[iii+1]]) * gca$full_moments[1:(iii+1)])) / gca$hn[iii+1] 
		retval <- retval + ci * gca$intpoly[[iii+1]](gca$x)
	}
	# sanity check; shall I throw a warning?
	retval <- pmin(1,pmax(0,retval))

	# support support
	if (is.finite(min(gca$support))) {
		retval[gca$x <= min(gca$support)] <- 0
	}
	if (is.finite(max(gca$support))) {
		retval[gca$x >= max(gca$support)] <- 1
	}

	# must be a better way to do these ... 
	if (!lower.tail) {
		retval <- 1 - retval
	}
	if (log.p)
		retval <- log(retval)
	return(retval)
}#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
