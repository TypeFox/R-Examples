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

# Created: 2015.03.16
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav
# Comments: Steven E. Pav

# transform the problem
dnf_tprob <- function(df1,df2,ncp1,ncp2) {
	n <- sum(mapply(function(d1,d2,n1,n2) { 1 },df1,df2,ncp1,ncp2,SIMPLIFY=TRUE))
	pow <- c(rep(1,n),rep(-1,n))
	df <- c(rep_len(df1,n),rep_len(df2,n))
	ncp <- c(rep_len(ncp1,n),rep_len(ncp2,n))
	retval <- list(n=n,pow=pow,df=df,ncp=ncp)
	return(retval)
}

# dproddnf, pproddnf, qproddnf, rproddnf#FOLDUP
#' @title The product of multiple doubly non-central F's distribution.
#'
#' @description 
#'
#' Density, distribution function, quantile function and random
#' generation for the product of multiple independent 
#' doubly non-central F variates.
#'
#' @details
#'
#' Let 
#' \deqn{x_j \sim F\left(\delta_{1,j},\delta_{2,j},\nu_{1,j},\nu_{2,j}\right)}{x_j ~ F(delta_1j,delta_2j,v_1j,v_2j)}
#' be independent doubly non-central F variates with non-centrality parameters
#' \eqn{\delta_{i,j}}{delta_ij} and degrees of freedom
#' \eqn{\nu_{i,j}}{v_ij} for \eqn{i=1,2,\ldots,I}{i=1,2,...,I} and 
#' \eqn{j=1,2}{j=1,2}.
#' Then 
#' \deqn{Y = \prod_j x_j}{Y = prod x_j}
#' takes a product of doubly non-central F's distribution. We take the
#' parameters of this distribution as the four \eqn{I}{I} length vectors
#' of the two degrees of freedom and the two non-centrality parameters.
#'
#' @usage
#'
#' dproddnf(x, df1, df2, ncp1, ncp2, log = FALSE, order.max=4)
#'
#' pproddnf(q, df1, df2, ncp1, ncp2, lower.tail = TRUE, log.p = FALSE, order.max=4)
#'
#' qproddnf(p, df1, df2, ncp1, ncp2, lower.tail = TRUE, log.p = FALSE, order.max=4)
#'
#' rproddnf(n, df1, df2, ncp1, ncp2)
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#'
#' @param df1,df2 the vectors of the degrees of freedom for the numerator and denominator.
#' We do \emph{not} recycle these versus the \code{x,q,p,n}.
#' @param ncp1,ncp2 the vectors of the non-centrality parameters for the numerator and denominator.
#' We do \emph{not} recycle these versus the \code{x,q,p,n}.
#'
#' @return \code{dproddnf} gives the density, \code{pproddnf} gives the 
#' distribution function, \code{qproddnf} gives the quantile function, 
#' and \code{rproddnf} generates random deviates.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @aliases dproddnf pproddnf qproddnf rproddnf
#' @seealso 
#' The sum of log of chi-squares distribution,
#' \code{\link{dsumlogchisq}},
#' \code{\link{psumlogchisq}},
#' \code{\link{qsumlogchisq}},
#' \code{\link{rsumlogchisq}}.
#' (doubly non-central) F distribution functions, 
#' \code{\link{ddnf}},
#' \code{\link{pdnf}}, 
#' \code{\link{qdnf}}, 
#' \code{\link{rdnf}}.
#'
#' @template etc
#' @template distribution
#' @template apx_distribution
#' @template not-recycled
#' @template ref-lnc
#' @note The PDQ functions are computed by translation of the 
#' sum of log chi-squares distribution functions.
#' @examples 
#' df1 <- c(10,20,5)
#' df2 <- c(1000,500,150)
#' ncp1 <- c(1,0,2.5)
#' ncp2 <- c(0,1.5,5)
#'
#' rv <- rproddnf(500, df1=df1,df2=df2,ncp1=ncp1,ncp2=ncp2)
#' d1 <- dproddnf(rv, df1=df1,df2=df2,ncp1=ncp1,ncp2=ncp2)
#' \dontrun{
#' plot(rv,d1)
#' }
#' p1 <- pproddnf(rv, df1=df1,df2=df2,ncp1=ncp1,ncp2=ncp2)
#' # should be nearly uniform:
#' \dontrun{
#' plot(ecdf(p1))
#' }
#' q1 <- qproddnf(ppoints(length(rv)), df1=df1,df2=df2,ncp1=ncp1,ncp2=ncp2)
#' \dontrun{
#' qqplot(x=rv,y=q1)
#' }
#' @name proddnf
#' @rdname dproddnf
#' @export 
dproddnf <- function(x, df1, df2, ncp1, ncp2, log = FALSE, order.max=4) {
	# beware bad recycling of the parameters!
	sub_p <- dnf_tprob(df1,df2,ncp1,ncp2)
	retval <- dsumlogchisq(log(x) + sum(sub_p$pow * log(sub_p$df)),
												 wts=sub_p$pow,df=sub_p$df,ncp=sub_p$ncp,log=log,order.max=order.max)
	if (log) {
		retval <- retval - log(x)
	} else {
		retval <- retval / x
	}
	return(retval)
}
#' @export 
pproddnf <- function(q, df1, df2, ncp1, ncp2, lower.tail = TRUE, log.p = FALSE, order.max=4) {
	# beware bad recycling of the parameters!
	sub_p <- dnf_tprob(df1,df2,ncp1,ncp2)
	retval <- psumlogchisq(log(q) + sum(sub_p$pow * log(sub_p$df)),
												 wts=sub_p$pow,df=sub_p$df,ncp=sub_p$ncp,
												 lower.tail=lower.tail,log.p=log.p,order.max=order.max)
	return(retval)
}
#' @export 
qproddnf <- function(p, df1, df2, ncp1, ncp2, lower.tail = TRUE, log.p = FALSE, order.max=4) {
	sub_p <- dnf_tprob(df1,df2,ncp1,ncp2)
	retval <- qsumlogchisq(p,wts=sub_p$pow,df=sub_p$df,ncp=sub_p$ncp,
												 lower.tail=lower.tail,log.p=log.p,order.max=order.max)
	retval <- exp(retval - sum(sub_p$pow * log(sub_p$df)))
	return(retval)
}
#' @export 
rproddnf <- function(n, df1, df2, ncp1, ncp2) {
	subX <- mapply(function(d1,d2,n1,n2)
										 { rdnf(n,df1=d1,df2=d2,ncp1=n1,ncp2=n2) },
										 df1,df2,ncp1,ncp2,SIMPLIFY=FALSE)
	X <- Reduce('*', subX)
	return(X)
}
#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
