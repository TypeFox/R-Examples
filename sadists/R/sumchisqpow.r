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

#
# Created: 2015.02.27
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav
# Comments: Steven E. Pav

# compute the cumulants of the sumchisqpow
# distribution. 
sumchisqpow_cumuls <- function(wts,df,ncp=0,pow=1,order.max=3) {
	#nterms <- max(vapply(list(wts,df,ncp),length,0))
	subkappa <- mapply(function(w,dd,nn,pp) 
										 { (w ^ (1:order.max)) * chipow_cumuls(df=dd,ncp=nn,pow=pp,order.max=order.max) },
										 wts,df,ncp,pow,SIMPLIFY=FALSE)
	kappa <- Reduce('+', subkappa)
	return(kappa)
}
sumchisqpow_support <- function(wts,df,ncp=0,pow=1) {
	minv <- ifelse(min(wts) < 0,-Inf,0)
	maxv <- ifelse(max(wts) > 0,Inf,0)
	retval <- c(minv,maxv)
	return(retval)
}

# dsumchisqpow, psumchisqpow, qsumchisqpow, rsumchisqpow#FOLDUP
#' @title The sum of (non-central) chi-squares raised to powers distribution.
#'
#' @description 
#'
#' Density, distribution function, quantile function and random
#' generation for the distribution of the weighted sum of non-central
#' chi-squares taken to powers.
#'
#' @details
#'
#' Let \eqn{X_i \sim \chi^2\left(\delta_i, \nu_i\right)}{X_i ~ chi^2(delta_i, v_i)}
#' be independently distributed non-central chi-squares, where \eqn{\nu_i}{v_i}
#' are the degrees of freedom, and \eqn{\delta_i}{delta_i} are the
#' non-centrality parameters.  
#' Let \eqn{w_i} and \eqn{p_i} be given constants. Suppose
#' \deqn{Y = \sum_i w_i X_i^{p_i}.}{Y = sum w_i (X_i)^(p_i).}
#' Then \eqn{Y}{Y} follows a weighted sum of chi-squares to power distribution. 
#'
#' @usage
#'
#' dsumchisqpow(x, wts, df, ncp=0, pow=1, log = FALSE, order.max=6)
#'
#' psumchisqpow(q, wts, df, ncp=0, pow=1, lower.tail = TRUE, log.p = FALSE, order.max=6)
#'
#' qsumchisqpow(p, wts, df, ncp=0, pow=1, lower.tail = TRUE, log.p = FALSE, order.max=6)
#'
#' rsumchisqpow(n, wts, df, ncp=0, pow=1)
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param wts the vector of weights. 
#' This is recycled against the \code{df, ncp, pow}, but not against the \code{x,q,p,n}.
#' @param df the vector of degrees of freedom. 
#' This is recycled against the \code{wts, ncp, pow}, but not against the \code{x,q,p,n}.
#' @param ncp the vector of non-centrality parameters. 
#' This is recycled against the \code{wts, df, pow}, but not against the \code{x,q,p,n}.
#' @param pow the vector of the power parameters. 
#' This is recycled against the \code{wts, df, ncp}, but not against the \code{x,q,p,n}.
#'
#' @template distribution
#' @template apx_distribution
#' @template not-recycled
#'
#' @return \code{dsumchisqpow} gives the density, \code{psumchisqpow} gives the 
#' distribution function, \code{qsumchisqpow} gives the quantile function, 
#' and \code{rsumchisqpow} generates random deviates.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @note The 'sum of chisquare power' distribution does \emph{not} generalize
#' the 'chi-bar-square' distribution, whose \emph{density} is the sum of
#' chi-square densities.
#' @aliases dsumchisqpow psumchisqpow qsumchisqpow rsumchisqpow
#' @seealso 
#' The upsilon distribution, 
#' \code{\link{dupsilon},\link{pupsilon},\link{qupsilon},\link{rupsilon}}.
#' @template etc
#' @examples 
#' wts <- c(1,-3,4)
#' df <- c(100,20,10)
#' ncp <- c(5,3,1)
#' pow <- c(1,0.5,1)
#' rvs <- rsumchisqpow(128, wts, df, ncp, pow)
#' dvs <- dsumchisqpow(rvs, wts, df, ncp, pow)
#' qvs <- psumchisqpow(rvs, wts, df, ncp, pow)
#' pvs <- qsumchisqpow(ppoints(length(rvs)), wts, df, ncp, pow)
#' @rdname dsumchisqpow
#' @name sumchisqpow
#' @export 
dsumchisqpow <- function(x, wts, df, ncp=0, pow=1, log = FALSE, order.max=6) {
	kappa <- sumchisqpow_cumuls(wts,df,ncp,pow,order.max=order.max)
	retval <- PDQutils::dapx_edgeworth(x,kappa,support=sumchisqpow_support(wts,df,ncp,pow),log=log)
	return(retval)
}
#' @export
psumchisqpow <- function(q, wts, df, ncp=0, pow=1, lower.tail = TRUE, log.p = FALSE, order.max=6) {
	kappa <- sumchisqpow_cumuls(wts,df,ncp,pow,order.max=order.max)
	retval <- PDQutils::papx_edgeworth(q,kappa,support=sumchisqpow_support(wts,df,ncp,pow),
																		 lower.tail=lower.tail,log.p=log.p)
	return(retval)
}
#' @export 
qsumchisqpow <- function(p, wts, df, ncp=0, pow=1, lower.tail = TRUE, log.p = FALSE, order.max=6) {
	kappa <- sumchisqpow_cumuls(wts,df,ncp,pow,order.max=order.max)
	retval <- PDQutils::qapx_cf(p,kappa,support=sumchisqpow_support(wts,df,ncp,pow),
															lower.tail=lower.tail,log.p=log.p)
	return(retval)
}
#' @export 
rsumchisqpow <- function(n, wts, df, ncp=0, pow=1) {
	subX <- mapply(function(w,dd,nn,pp) { w * (unbroken_rchisq(n,df=dd,ncp=nn) ^ pp) },
										 wts,df,ncp,pow,SIMPLIFY=FALSE)
	X <- Reduce('+', subX)
	return(X)
}
#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
