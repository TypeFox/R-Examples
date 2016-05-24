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

#source("utils.r")

# compute the moments of the doubly non-central
# t distribution. this is distributed as
#
# (ncp_1 + Z) / sqrt(chi^2_2(ncp_2) / df_2)
#
# where the Z is independent of the chi-square.
dnt_moments <- function(df,ncp1,ncp2,order.max=3) {
	orders = 1:order.max
	numr <- norm_moms(mu=ncp1,sigma=1,order.max=order.max)
	deno <- exp(chisq_moms(df=df,ncp=ncp2,orders=-orders/2.0,log=TRUE) + (orders/2.0)*log(df))
	mus <- numr * deno
	return(mus)
}

# compute the moments of the doubly non-central
# t distribution. this is distributed as
#
# (ncp_1 + Z) / sqrt(chi^2_2(ncp_2) / df_2)
#
# where the Z is independent of the chi-square.
dnt_cumuls <- function(df,ncp1,ncp2,order.max=3) {
	kappa <- moment2cumulant(dnt_moments(df,ncp1,ncp2,order.max))
	return(kappa)
}

# ddnt, pdnt, qdnt, rdnt#FOLDUP
#' @title The doubly non-central t distribution.
#'
#' @description 
#'
#' Density, distribution function, quantile function and random
#' generation for the doubly non-central t distribution.
#'
#' @details
#'
#' Let \eqn{Z \sim \mathcal{N}\left(\mu,1\right)}{Z ~ N(u,1)} independently
#' of \eqn{X \sim \chi^2\left(\theta,\nu\right)}{X ~ x^2(theta,v)}. The 
#' random variable
#' \deqn{T = \frac{Z}{\sqrt{X/\nu}}}{T = Z / sqrt(X/v)}
#' takes a \emph{doubly non-central t distribution} with parameters
#' \eqn{\nu, \mu, \theta}{v, mu, theta}.
#'
#' @usage
#'
#' ddnt(x, df, ncp1, ncp2, log = FALSE, order.max=6)
#'
#' pdnt(q, df, ncp1, ncp2, lower.tail = TRUE, log.p = FALSE, order.max=6)
#'
#' qdnt(p, df, ncp1, ncp2, lower.tail = TRUE, log.p = FALSE, order.max=6)
#'
#' rdnt(n, df, ncp1, ncp2)
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#'
#' @param df the degrees of freedom for the denominator, \eqn{\nu}{v}.
#' We do \emph{not} recycle these versus the \code{x,q,p,n}.
#' @param ncp1,ncp2 the non-centrality parameters for the numerator and denominator,
#' respectively, \eqn{\mu}{mu} and \eqn{\theta}{theta}
#' We do \emph{not} recycle these versus the \code{x,q,p,n}.
#' 
#' @return \code{ddnt} gives the density, \code{pdnt} gives the 
#' distribution function, \code{qdnt} gives the quantile function, 
#' and \code{rdnt} generates random deviates.
#' 
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @aliases ddnt pdnt qdnt rdnt
#' @seealso t distribution functions, \code{\link{dt}, \link{pt}, \link{qt}, \link{rt}}
#' @template etc
#' @template distribution
#' @template apx_distribution
#' @template not-recycled
#' @template dnt
#' @name dnt 
#' @rdname ddnt
#' @examples 
#' rvs <- rdnt(128, 20, 1, 1)
#' dvs <- ddnt(rvs, 20, 1, 1)
#' pvs.H0 <- pdnt(rvs, 20, 0, 1)
#' pvs.HA <- pdnt(rvs, 20, 1, 1)
#' \dontrun{
#' plot(ecdf(pvs.H0))
#' plot(ecdf(pvs.HA))
#' }
#' # compare to singly non-central
#' dv1 <- ddnt(1, df=10, ncp1=5, ncp2=0, log=FALSE)
#' dv2 <- dt(1, df=10, ncp=5, log=FALSE)
#' pv1 <- pdnt(1, df=10, ncp1=5, ncp2=0, log.p=FALSE)
#' pv11 <- pdnt(1, df=10, ncp1=5, ncp2=0.001, log.p=FALSE)
#' v2 <- pt(1, df=10, ncp=5, log.p=FALSE)
#'
#' q1 <- qdnt(pv1, df=10, ncp1=5, ncp2=0, log.p=FALSE)
#' @export
ddnt <- function(x,df,ncp1,ncp2,log=FALSE, order.max=6) {
	kappa <- dnt_cumuls(df,ncp1,ncp2,order.max=order.max)
	retval <- PDQutils::dapx_edgeworth(x,kappa,log=log)
	return(retval)
}
#' @export
pdnt <- function(q, df, ncp1, ncp2, lower.tail = TRUE, log.p = FALSE, order.max=6) {
	kappa <- dnt_cumuls(df,ncp1,ncp2,order.max=order.max)
	retval <- PDQutils::papx_edgeworth(q,kappa,lower.tail=lower.tail,log.p=log.p)
	return(retval)
}
#' @export
qdnt <- function(p, df, ncp1, ncp2, lower.tail = TRUE, log.p = FALSE, order.max=6) {
	kappa <- dnt_cumuls(df,ncp1,ncp2,order.max=order.max)
	retval <- PDQutils::qapx_cf(p,kappa,lower.tail=lower.tail,log.p=log.p)
	return(retval)
}
#' @export 
rdnt <- function(n,df,ncp1,ncp2) {
	X <- rnorm(n,mean=ncp1,sd=1)
	Y <- rchisq(n,df=df,ncp=ncp2)
	Z <- X / sqrt(Y / df)
	return(Z)
}
#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
