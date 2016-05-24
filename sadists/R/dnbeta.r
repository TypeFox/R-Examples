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

# ddnbeta, pdnbeta, qdnbeta, rdnbeta#FOLDUP
#' @title The doubly non-central Beta distribution.
#'
#' @description 
#'
#' Density, distribution function, quantile function and random
#' generation for the doubly non-central Beta distribution.
#'
#' @details
#'
#' Suppose \eqn{x_i \sim \chi^2\left(\delta_i,\nu_i\right)}{x_i ~ X^2(delta_i,v_i)}
#' be independent non-central chi-squares for \eqn{i=1,2}{i=1,2}.
#' Then 
#' \deqn{Y = \frac{x_1}{x_1 + x_2}}{Y = x_1 / (x_1 + x_2)}
#' takes a doubly non-central Beta distribution with degrees of freedom
#' \eqn{\nu_1, \nu_2}{v_1, v_2} and non-centrality parameters
#' \eqn{\delta_1,\delta_2}{delta_1,delta_2}.
#'
#' @usage
#'
#' ddnbeta(x, df1, df2, ncp1, ncp2, log = FALSE, order.max=6)
#'
#' pdnbeta(q, df1, df2, ncp1, ncp2, lower.tail = TRUE, log.p = FALSE, order.max=6)
#'
#' qdnbeta(p, df1, df2, ncp1, ncp2, lower.tail = TRUE, log.p = FALSE, order.max=6)
#'
#' rdnbeta(n, df1, df2, ncp1, ncp2)
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#'
#' @param df1,df2 the degrees of freedom for the numerator and denominator.
#' We do \emph{not} recycle these versus the \code{x,q,p,n}.
#' @param ncp1,ncp2 the non-centrality parameters for the numerator and denominator.
#' We do \emph{not} recycle these versus the \code{x,q,p,n}.
#'
#' @return \code{ddnbeta} gives the density, \code{pdnbeta} gives the 
#' distribution function, \code{qdnbeta} gives the quantile function, 
#' and \code{rdnbeta} generates random deviates.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @aliases ddnbeta pdnbeta qdnbeta rdnbeta
#' @seealso (doubly non-central) F distribution functions, 
#' \code{\link{ddnf}, \link{pdnf}, \link{qdnf}, \link{rdnf}}.
#' @template etc
#' @template distribution
#' @template apx_distribution
#' @template not-recycled
#' @examples 
#' rv <- rdnbeta(500, df1=100,df2=500,ncp1=1.5,ncp2=12)
#' d1 <- ddnbeta(rv, df1=100,df2=500,ncp1=1.5,ncp2=12)
#' \dontrun{
#' plot(rv,d1)
#' }
#' p1 <- ddnbeta(rv, df1=100,df2=500,ncp1=1.5,ncp2=12)
#' # should be nearly uniform:
#' \dontrun{
#' plot(ecdf(p1))
#' }
#' q1 <- qdnbeta(ppoints(length(rv)), df1=100,df2=500,ncp1=1.5,ncp2=12)
#' \dontrun{
#' qqplot(x=rv,y=q1)
#' }
#' @name dnbeta
#' @rdname ddnbeta
#' @export 
ddnbeta <- function(x, df1, df2, ncp1, ncp2, log = FALSE, order.max=6) {
	xF <- (df2/df1) * x / (1-x)
	retval <- ddnf(xF,df1=df1,df2=df2,ncp1=ncp1,ncp2=ncp2,log=log,order.max=order.max)
	if (log) {
		retval <- log(df2/df1) + retval - 2 * log(1-x)
	} else {
		retval <- (df2/df1) * retval / ((1-x)^2)
	}
	in_range <- 0 <= x & x < 1
	retval[!in_range] <- NaN
	return(retval)
}
#' @export 
pdnbeta <- function(q, df1, df2, ncp1, ncp2, lower.tail = TRUE, log.p = FALSE, order.max=6) {
	qF <- (df2/df1) * q / (1-q)
	retval <- pdnf(qF,df1=df1,df2=df2,ncp1=ncp1,ncp2=ncp2,lower.tail=lower.tail,log.p=log.p,order.max=order.max)
	minv <- 1 - as.double(lower.tail)
	maxv <- 1 - minv
	if (log.p) {
		minv <- log(minv)
		maxv <- log(maxv)
	}
	retval[q < 0] <- minv
	retval[q >= 1] <- maxv
	return(retval)
}
#' @export 
qdnbeta <- function(p, df1, df2, ncp1, ncp2, lower.tail = TRUE, log.p = FALSE, order.max=6) {
	qF <- (df1/df2) * qdnf(p,df1=df1,df2=df2,ncp1=ncp1,ncp2=ncp2,lower.tail=lower.tail,log.p=log.p,order.max=order.max)
	retval <- qF / (1 + qF)
	in_range <- ifelse(log.p,p <= 0,0 <= p & p <= 1)
	retval[!in_range] <- NaN
	return(retval)
}
#' @export 
rdnbeta <- function(n, df1, df2, ncp1, ncp2) {
	X1 <- unbroken_rchisq(n,df=df1,ncp=ncp1)
	X2 <- unbroken_rchisq(n,df=df2,ncp=ncp2)
	X <- X1 / (X1+X2)
	return(X)
}
#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
