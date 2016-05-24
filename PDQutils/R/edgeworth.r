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

# Created: 2015.02.19
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav
# Comments: Steven E. Pav

# for the Hermite Polynomials
require(orthopolynom)
require(moments)

# blinnikov and moessner's ADVANCE function
advance <- function(kms) {
	current <- kms$current
	mold <- kms$mold

	n <- length(current)
	ords <- 1:(length(current))
	stopifnot(n == sum(current * ords))
	sumcur <- n
	m <- 1
	is.done <- FALSE
	while (!is.done) {
		sumcur <- sumcur - current[m]*m + (m+1)
		current[m] <- 0
		current[m+1] <- current[m+1] + 1
		m <- m+1
		is.done <- (sumcur <= n) || (m > mold)
	}
	mold <- max(mold,m)
	current[1] <- n - sumcur
	retv <- list(current=current,mold=mold)
	return(retv)
}

#' @title Approximate density and distribution via Edgeworth expansion.
#'
#' @description 
#'
#' Approximate the probability density or cumulative distribution function of a distribution via its raw cumulants.
#'
#' @template details-edgeworth 
#'
#' @usage
#'
#' dapx_edgeworth(x, raw.cumulants, support=c(-Inf,Inf), log=FALSE)
#'
#' papx_edgeworth(q, raw.cumulants, support=c(-Inf,Inf), lower.tail=TRUE, log.p=FALSE)
#'
#' @param x where to evaluate the approximate density.
#' @param q where to evaluate the approximate distribution.
#' @param raw.cumulants an atomic array of the 1st through kth raw cumulants
#' of the probability distribution. The first cumulant is the mean, the
#' second is the variance. The third is \emph{not} the typical unitless skew.
#' @param support the support of the density function. It is assumed
#' that the density is zero on the complement of this open interval.
#' @param log logical; if TRUE, densities \eqn{f} are given 
#'  as \eqn{\mbox{log}(f)}{log(f)}.
#' @param log.p logical; if TRUE, probabilities p are given 
#'  as \eqn{\mbox{log}(p)}{log(p)}.
#' @param lower.tail whether to compute the lower tail. If false, we approximate the survival function.
#' @return The approximate density at \code{x}, or the approximate CDF at
#' \code{q}.
#'
#' @keywords distribution 
#' @seealso the Gram Charlier expansions, \code{\link{dapx_gca}, \link{papx_gca}}
#' @export 
#' @template ref-Blinnikov
#' @aliases papx_edgeworth 
#' @note 
#'
#' Monotonicity of the CDF is not guaranteed.
#'
#' @examples 
#' # normal distribution, for which this is silly
#' xvals <- seq(-2,2,length.out=501)
#' d1 <- dapx_edgeworth(xvals, c(0,1,0,0,0,0))
#' d2 <- dnorm(xvals)
#' d1 - d2
#'
#' qvals <- seq(-2,2,length.out=501)
#' p1 <- papx_edgeworth(qvals, c(0,1,0,0,0,0))
#' p2 <- pnorm(qvals)
#' p1 - p2
#' @template etc
dapx_edgeworth <- function(x,raw.cumulants,support=c(-Inf,Inf),log=FALSE) {#FOLDUP
	order.max <- length(raw.cumulants)

	mu <- raw.cumulants[1]
	sigma <- sqrt(raw.cumulants[2])

	eta <- (x - mu) / sigma

	phi.eta <- dnorm(eta)
	retval <- phi.eta

	# compute via equation (43) of Blinnikov and Moessner
	if (order.max > 2) {
		hermi <- orthopolynom::hermite.he.polynomials(3*order.max+2, normalized=FALSE)
		Sn <- raw.cumulants[3:order.max] / (raw.cumulants[2] ^ (2:(order.max-1)))

		for (s in c(1:(order.max-2))) {
			nexterm <- rep(0,length(phi.eta))
			kms <- list(mold=1,current=c(s,rep(0,s-1)))
			while (kms$mold <= s) {
				r <- sum(kms$current)
				coefs <- ((Sn[1:s] / factorial(3:(s+2))) ^ kms$current) / factorial(kms$current)
				coef <- prod(coefs)
				nexterm <- nexterm + coef * as.function(hermi[[s+2*r+1]])(eta)
				kms <- advance(kms)
			}
			retval <- retval + phi.eta * (sigma^s) * nexterm
		}
	}

	# adjust back from standardized
	retval <- retval / sigma

	# sanity check; shall I throw a warning?
	retval <- pmax(0,retval)

	# support support
	if (is.finite(min(support))) {
		retval[x <= min(support)] <- 0
	}
	if (is.finite(max(support))) {
		retval[x >= max(support)] <- 0
	}

	# must be a better way to do this ... 
	if (log)
		retval <- log(retval)
	return(retval)
}#UNFOLD
#' @export
papx_edgeworth <- function(q,raw.cumulants,support=c(-Inf,Inf),lower.tail=TRUE,log.p=FALSE) {#FOLDUP
	order.max <- length(raw.cumulants)

	# 2FIX: would it not be better to pass lower.tail to pnorm and dnorm below?
	# or subtract 1 from the end result?
	if (!lower.tail) {
		# transform q and the raw cumulants
		q <- - q
		raw.cumulants <- raw.cumulants * ((-1)^(1:order.max))
		support <- sort(-support)
	}

	mu <- raw.cumulants[1]
	sigma <- sqrt(raw.cumulants[2])

	eta <- (q - mu) / sigma

	phi.eta <- dnorm(eta)
	retval <- pnorm(eta)

	# compute via equation (43) of Blinnikov and Moessner
	if (order.max > 2) {
		hermi <- orthopolynom::hermite.he.polynomials(3*order.max+2, normalized=FALSE)
		Sn <- raw.cumulants[3:order.max] / (raw.cumulants[2] ^ (2:(order.max-1)))

		for (s in c(1:(order.max-2))) {
			nexterm <- rep(0,length(phi.eta))
			kms <- list(mold=1,current=c(s,rep(0,s-1)))
			while (kms$mold <= s) {
				r <- sum(kms$current)
				coefs <- ((Sn[1:s] / factorial(3:(s+2))) ^ kms$current) / factorial(kms$current)
				coef <- prod(coefs)
				# n.b. the hermite polynomial is one order less than in the dapx
				# and we _subtract_ it
				nexterm <- nexterm - coef * as.function(hermi[[s+2*r]])(eta)
				kms <- advance(kms)
			}
			retval <- retval + phi.eta * (sigma^s) * nexterm
		}
	}

	# sanity check; shall I throw a warning?
	retval <- pmin(1,pmax(0,retval))

	# support support
	if (is.finite(min(support))) {
		retval[q <= min(support)] <- 0
	}
	if (is.finite(max(support))) {
		retval[q >= max(support)] <- 1
	}

	# must be a better way to do this ... 
	if (log.p)
		retval <- log(retval)
	return(retval)
}#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
