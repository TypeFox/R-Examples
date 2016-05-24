# Copyright 2014-2014 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav
#
# This file is part of sadists.
#
# sadists is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# sadists is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with sadists.  If not, see <http://www.gnu.org/licenses/>.

# Created: 2014.02.14
# Copyright: Steven E. Pav, 2014
# Author: Steven E. Pav
# Comments: Steven E. Pav

# compute the raw moments of the K-prime
# distribution. this is distributed as
# (b Z + a sqrt(chi^2_v1 / v1)) / (sqrt(chi^2_v2 /v2))
# or 
# b * (Z + (a/b) sqrt(chi^2_v1/v1)) / (sqrt(chi^2_v2 / v2))
# or
# (b * sqrt(v2)) * lambda_prime(v1,a/b) / (sqrt(chi^2_v2))
#
# when b == 0 this is
# a sqrt(v2/v1) * sqrt(chi^2_v1 / chi^2_v2)
kprime_moms <- function(v1,v2,a,b,order.max=3) {
	orders <- 1:order.max
	if (b != 0) {
		mu <- ((b * sqrt(v2)) ^ (orders)) * lambdap_moms(df=v1,t=a/b,order.max=order.max) * 
						chisq_moms(df=v2,ncp=0,orders=-orders/2.0)
	} else {
		# 2FIX: use logs for stability?
		mu <- ((a*sqrt(v2/v1))^orders) * chisq_moms(df=v1,ncp=0,orders=orders/2.0) *
						chisq_moms(df=v2,ncp=0,orders=-orders/2.0)
	}
	return(mu)
}
# compute the raw cumulants of the K-prime
# distribution. this is distributed as
# (b Z + a sqrt(chi^2_v1 / v1)) / (sqrt(chi^2_v2 /v2))
kprime_cumuls <- function(v1,v2,a,b,order.max=3) {
	kappa <- moment2cumulant(kprime_moms(v1,v2,a,b,order.max))
	return(kappa)
}

# dkprime, pkprime, qkprime, rkprime#FOLDUP
#' @title The K prime distribution.
#'
#' @description 
#'
#' Density, distribution function, quantile function and random
#' generation for the K prime distribution.
#'
#' @details
#'
#' Suppose \eqn{y \sim \chi^2\left(\nu_1\right)}{y ~ x^2(v1)}, and
#' \eqn{x \sim t \left(\nu_2, a\sqrt{y/\nu_1}/b\right)}{x ~ t(v2,(a/b) sqrt(y/v1))}.
#' Then the random variable
#' \deqn{T = b x}{T = b x}
#' takes a K prime distribution with parameters 
#' \eqn{\nu_1, \nu_2, a, b}{v1, v2, a, b}. In Lecoutre's terminology,
#' \eqn{T \sim K'_{\nu_1, \nu_2}\left(a, b\right)}{T ~ K'_v1,v2(a,b)}
#'
#' Equivalently, we can think of
#' \deqn{T = \frac{b Z + a \sqrt{\chi^2_{\nu_1} / \nu_1}}{\sqrt{\chi^2_{\nu_2} / \nu_2}}}{T = (bZ + a sqrt(chi2_v1/v1)) / sqrt(chi2_v2/v2)}
#' where \eqn{Z} is a standard normal, and the normal and the (central) chi-squares are
#' independent of each other. When \eqn{a=0}{a=0} we recover
#' a central t distribution; 
#' when \eqn{\nu_1=\infty}{v1=inf} we recover a rescaled non-central t distribution;
#' when \eqn{b=0}{b=0}, we get a rescaled square root of a central F
#' distribution; when \eqn{\nu_2=\infty}{v2=inf}, we recover a 
#' Lambda prime distribution.
#'
#' @usage
#'
#' dkprime(x, v1, v2, a, b = 1, order.max=6, log = FALSE)
#'
#' pkprime(q, v1, v2, a, b = 1, order.max=6, lower.tail = TRUE, log.p = FALSE)
#'
#' qkprime(p, v1, v2, a, b = 1, order.max=6, lower.tail = TRUE, log.p = FALSE)
#'
#' rkprime(n, v1, v2, a, b = 1)
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param v1 the degrees of freedom in the numerator chisquare. When
#' (positive) infinite, we recover a non-central t 
#' distribution with \code{v2} degrees of freedom and non-centrality
#' parameter \code{a}, scaled by \code{b}.
#' This is not recycled against the \code{x,q,p,n}.
#' @param v2 the degrees of freedom in the denominator chisquare.
#' When equal to infinity, we recover the Lambda prime distribution.
#' This is not recycled against the \code{x,q,p,n}.
#' @param a the non-centrality scaling parameter. When equal to zero,
#' we recover the (central) t distribution.
#' This is not recycled against the \code{x,q,p,n}.
#' @param b the scaling parameter.
#' This is not recycled against the \code{x,q,p,n}.
#'
#' @return \code{dkprime} gives the density, \code{pkprime} gives the 
#' distribution function, \code{qkprime} gives the quantile function, 
#' and \code{rkprime} generates random deviates.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @aliases dkprime pkprime qkprime rkprime
#' @seealso t distribution functions, \code{\link{dt}, \link{pt}, \link{qt}, \link{rt}},
#' lambda prime distribution functions, \code{\link{dlambdap}, \link{plambdap}, \link{qlambdap}, \link{rlambdap}}.
#' @template distribution
#' @template apx_distribution
#' @template not-recycled
#' @template etc
#' @template ref-kprime
#' @examples 
#' d1 <- dkprime(1, 50, 20, a=0.01)
#' d2 <- dkprime(1, 50, 20, a=0.0001)
#' d3 <- dkprime(1, 50, 20, a=0)
#' d4 <- dkprime(1, 10000, 20, a=1)
#' d5 <- dkprime(1, Inf, 20, a=1)
#'
#' @rdname dkprime
#' @name kprime
#' @export 
dkprime <- function(x, v1, v2, a, b=1, order.max=6, log = FALSE) {
	kappa <- kprime_cumuls(v1,v2,a,b,order.max=order.max)
	retval <- PDQutils::dapx_edgeworth(x,kappa,log=log)
	return(retval)
}
#' @export 
pkprime <- function(q, v1, v2, a, b=1, order.max=6, lower.tail = TRUE, log.p = FALSE) {
	kappa <- kprime_cumuls(v1,v2,a,b,order.max=order.max)
	retval <- PDQutils::papx_edgeworth(q,kappa,lower.tail=lower.tail,log.p=log.p)
	return(retval)
}
#' @export 
qkprime <- function(p, v1, v2, a, b=1, order.max=6, lower.tail = TRUE, log.p = FALSE) {
	kappa <- kprime_cumuls(v1,v2,a,b,order.max=order.max)
	retval <- PDQutils::qapx_cf(p,kappa,lower.tail=lower.tail,log.p=log.p)
	return(retval)
}
#' @export 
rkprime <- function(n, v1, v2, a, b = 1) {
	#2FIX: check for b = 0...
	#y <- rchisq(n,df=v1) 
	#ncp <- sqrt(y/v1) * (a/b)
	#X <- b * rt(n,df=v2,ncp=ncp)
	X <- (rnorm(n,mean=0,sd=b) + a * sqrt(rchisq(n,df=v1)/v1)) / sqrt(rchisq(n,df=v2)/v2)
	return(X)
}
#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
