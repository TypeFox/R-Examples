# Copyright 2014-2015 Steven E. Pav. All Rights Reserved.
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

# Created: 2015.03.07
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav
# Comments: Steven E. Pav

# ddneta, pdneta, qdneta, rdneta#FOLDUP
#' @title The doubly non-central Eta distribution.
#'
#' @description 
#'
#' Density, distribution function, quantile function and random
#' generation for the doubly non-central Eta distribution.
#'
#' @details
#'
#' Suppose \eqn{Z}{Z} is a normal with mean \eqn{\delta_1}{delta_1},
#' independent of \eqn{X \sim \chi^2\left(\delta_2,\nu_2\right)}{X ~ X^2(delta_2,v_2)},
#' a non-central chi-square with \eqn{\nu_2}{v_2} degrees of freedom
#' and non-centrality parameter \eqn{\delta_2}{delta_2}. Then
#' \deqn{Y = \frac{Z}{\sqrt{Z^2 + X}}}{Y = Z/sqrt(Z^2 + X)}
#' takes a doubly non-central Eta distribution with 
#' \eqn{\nu_2}{v_2} degrees of freedom and non-centrality parameters
#' \eqn{\delta_1,\delta_2}{delta_1,delta_2}. The \emph{square} of
#' a doubly non-central Eta is a doubly non-central Beta variate.
#'
#' @usage
#'
#' ddneta(x, df, ncp1, ncp2, log = FALSE, order.max=6)
#'
#' pdneta(q, df, ncp1, ncp2, lower.tail = TRUE, log.p = FALSE, order.max=6)
#'
#' qdneta(p, df, ncp1, ncp2, lower.tail = TRUE, log.p = FALSE, order.max=6)
#'
#' rdneta(n, df, ncp1, ncp2)
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#'
#' @param df the degrees of freedom for the denominator chi square.
#' We do \emph{not} recycle this versus the \code{x,q,p,n}.
#' @param ncp1,ncp2 the non-centrality parameters for the numerator and denominator.
#' We do \emph{not} recycle these versus the \code{x,q,p,n}.
#'
#' @return \code{ddneta} gives the density, \code{pdneta} gives the 
#' distribution function, \code{qdneta} gives the quantile function, 
#' and \code{rdneta} generates random deviates.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @aliases ddneta pdneta qdneta rdneta
#' @seealso (doubly non-central) t distribution functions, 
#' \code{\link{ddnt}, \link{pdnt}, \link{qdnt}, \link{rdnt}}.
#' @seealso (doubly non-central) Beta distribution functions, 
#' \code{\link{ddnbeta}, \link{pdnbeta}, \link{qdnbeta}, \link{rdnbeta}}.
#' @template etc
#' @template distribution
#' @template apx_distribution
#' @template not-recycled
#' @examples 
#' rv <- rdneta(500, df=100,ncp1=1.5,ncp2=12)
#' d1 <- ddneta(rv, df=100,ncp1=1.5,ncp2=12)
#' \dontrun{
#' plot(rv,d1)
#' }
#' p1 <- ddneta(rv, df=100,ncp1=1.5,ncp2=12)
#' # should be nearly uniform:
#' \dontrun{
#' plot(ecdf(p1))
#' }
#' q1 <- qdneta(ppoints(length(rv)), df=100,ncp1=1.5,ncp2=12)
#' \dontrun{
#' qqplot(x=rv,y=q1)
#' }
#' @name dneta
#' @rdname ddneta
#' @export 
ddneta <- function(x, df, ncp1, ncp2, log = FALSE, order.max=6) {
	xF <- sqrt(df) * x / sqrt(1-x^2)
	retval <- ddnt(xF,df=df,ncp1=ncp1,ncp2=ncp2,log=log,order.max=order.max)
	if (log) {
		retval <- 0.5 * log(df) + retval - (1.5) * log(1-x^2)
	} else {
		retval <- sqrt(df) * retval / ((1-x^2)^(1.5))
	}
	in_range <- -1 <= x & x < 1
	retval[!in_range] <- NaN
	return(retval)
}
#' @export 
pdneta <- function(q, df, ncp1, ncp2, lower.tail = TRUE, log.p = FALSE, order.max=6) {
	qF <- sqrt(df) * q / sqrt(1-q^2)
	retval <- pdnt(qF,df=df,ncp1=ncp1,ncp2=ncp2,lower.tail=lower.tail,log.p=log.p,order.max=order.max)
	minv <- 1 - as.double(lower.tail)
	maxv <- 1 - minv
	if (log.p) {
		minv <- log(minv)
		maxv <- log(maxv)
	}
	retval[q < -1] <- minv
	retval[q >= 1] <- maxv
	return(retval)
}
#' @export 
qdneta <- function(p, df, ncp1, ncp2, lower.tail = TRUE, log.p = FALSE, order.max=6) {
# 2FIX: this part
	qF <- (1/sqrt(df)) * qdnt(p,df=df,ncp1=ncp1,ncp2=ncp2,lower.tail=lower.tail,log.p=log.p,order.max=order.max)
	retval <- qF / sqrt(1 + qF^2)
	in_range <- ifelse(log.p,p <= 0,0 <= p & p <= 1)
	retval[!in_range] <- NaN
	return(retval)
}
#' @export 
rdneta <- function(n, df, ncp1, ncp2) {
	X1 <- rnorm(n,mean=ncp1)
	X2 <- unbroken_rchisq(n,df=df,ncp=ncp2)
	X <- X1 / sqrt((X1^2) + X2)
	return(X)
}
#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
