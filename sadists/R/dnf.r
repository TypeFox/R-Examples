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

# Created: 2015.02.14
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav
# Comments: Steven E. Pav

# see also:
#
# http://www.statsresearch.co.nz/robert/QF.htm
# http://www.jstor.org/stable/2347725
# Pan's algo

# compute the moments of the doubly non-central
# F distribution. this is distributed as
#
# (chi^2_1(ncp_1) / df_1) / (chi^2_2(ncp_2) / df_2)
#
# where the chi^2 are independent chi-square
dnf_moments <- function(df1,df2,ncp1,ncp2,order.max=3) {
	orders = 1:order.max
	log.numr <- chisq_moms(df=df1,ncp=ncp1,orders=orders,log=TRUE) - orders*log(df1)
	log.deno <- chisq_moms(df=df2,ncp=ncp2,orders=-orders,log=TRUE) + orders*log(df2)
	mus <- exp(log.numr + log.deno)
	return(mus)
}

# compute the cumulants of the doubly non-central
# F distribution. this is distributed as
#
# (chi^2_1(ncp_1) / df_1) / (chi^2_2(ncp_2) / df_2)
#
# where the chi^2 are independent chi-square
dnf_cumuls <- function(df1,df2,ncp1,ncp2,order.max=3) {
	kappa <- moment2cumulant(dnf_moments(df1,df2,ncp1,ncp2,order.max))
	return(kappa)
}

# ddnf, pdnf, qdnf, rdnf#FOLDUP
#' @title The doubly non-central F distribution.
#'
#' @description 
#'
#' Density, distribution function, quantile function and random
#' generation for the doubly non-central F distribution.
#'
#' @details
#'
#' Suppose \eqn{x_i \sim \chi^2\left(\delta_i,\nu_i\right)}{x_i ~ X^2(delta_i,v_i)}
#' be independent non-central chi-squares for \eqn{i=1,2}{i=1,2}.
#' Then 
#' \deqn{Y = \frac{x_1/\nu_1}{x_2/\nu_2}}{Y = (x_1/v_1) / (x_2/v_2)}
#' takes a doubly non-central F distribution with degrees of freedom
#' \eqn{\nu_1, \nu_2}{v_1, v_2} and non-centrality parameters
#' \eqn{\delta_1,\delta_2}{delta_1,delta_2}.
#'
#' @usage
#'
#' ddnf(x, df1, df2, ncp1, ncp2, log = FALSE, order.max=6)
#'
#' pdnf(q, df1, df2, ncp1, ncp2, lower.tail = TRUE, log.p = FALSE, order.max=6)
#'
#' qdnf(p, df1, df2, ncp1, ncp2, lower.tail = TRUE, log.p = FALSE, order.max=6)
#'
#' rdnf(n, df1, df2, ncp1, ncp2)
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
#' @return \code{ddnf} gives the density, \code{pdnf} gives the 
#' distribution function, \code{qdnf} gives the quantile function, 
#' and \code{rdnf} generates random deviates.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @aliases ddnf pdnf qdnf rdnf
#' @seealso (singly non-central) F distribution functions, 
#' \code{\link{df}, \link{pf}, \link{qf}, \link{rf}}.
#' @template etc
#' @template distribution
#' @template apx_distribution
#' @template not-recycled
#' @examples 
#' rv <- rdnf(500, df1=100,df2=500,ncp1=1.5,ncp2=12)
#' d1 <- ddnf(rv, df1=100,df2=500,ncp1=1.5,ncp2=12)
#' \dontrun{
#' plot(rv,d1)
#' }
#' p1 <- ddnf(rv, df1=100,df2=500,ncp1=1.5,ncp2=12)
#' # should be nearly uniform:
#' \dontrun{
#' plot(ecdf(p1))
#' }
#' q1 <- qdnf(ppoints(length(rv)), df1=100,df2=500,ncp1=1.5,ncp2=12)
#' \dontrun{
#' qqplot(x=rv,y=q1)
#' }
#' @name dnf
#' @rdname ddnf
#' @export 
ddnf <- function(x, df1, df2, ncp1, ncp2, log = FALSE, order.max=6) {
	kappa <- dnf_cumuls(df1,df2,ncp1,ncp2,order.max=order.max)
	retval <- PDQutils::dapx_edgeworth(x,kappa,support=c(0,Inf),log=log)
	return(retval)
}
#' @export 
pdnf <- function(q, df1, df2, ncp1, ncp2, lower.tail = TRUE, log.p = FALSE, order.max=6) {
	kappa <- dnf_cumuls(df1,df2,ncp1,ncp2,order.max=order.max)
	retval <- PDQutils::papx_edgeworth(q,kappa,support=c(0,Inf),lower.tail=lower.tail,log.p=log.p)
	return(retval)
}
#' @export 
qdnf <- function(p, df1, df2, ncp1, ncp2, lower.tail = TRUE, log.p = FALSE, order.max=6) {
	kappa <- dnf_cumuls(df1,df2,ncp1,ncp2,order.max=order.max)
	retval <- PDQutils::qapx_cf(p,kappa,support=c(0,Inf),lower.tail=lower.tail,log.p=log.p)
	return(retval)
}
#' @export 
rdnf <- function(n, df1, df2, ncp1, ncp2) {
	X <- (rchisq(n,df=df1,ncp=ncp1)/df1) / (rchisq(n,df=df2,ncp=ncp2)/df2)
	return(X)
}
#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
