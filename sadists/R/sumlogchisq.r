# /usr/bin/r
#
# Created: 2015.03.18
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav <steven@cerebellumcapital.com>
# Comments: Steven E. Pav

# cumulants of the log of the central chi-square
lc_cumuls <- function(df,order.max=3,orders=c(1:order.max)) {
	kappa <- psigamma(df/2,deriv=orders-1)
	kappa[1] <- kappa[1] + log(2)
	return(kappa)
}
# moments of the log of the central chi-square
lc_moments <- function(df,order.max=3,orders=c(1:order.max)) {
	kappa <- lc_cumuls(df,orders=orders)
	mu <- PDQutils::cumulant2moment(kappa)
	return(mu)
}
# moments of the log of the non-central chi-square
lnc_moments <- function(df,ncp=0,order.max=3,orders=c(1:order.max)) {
	stopifnot(ncp >=0)
	if (ncp > 0) {
		hancp <- ncp / 2.0
		# should be smarter about 0:50 here.
		allmu <- sapply(0:50,function(iv) {
										exp(-hancp + iv * log(hancp) - lfactorial(iv)) * lc_moments(df+2*iv,orders=orders)
			},simplify=FALSE)
		mu <- Reduce('+',allmu)
	} else {
		mu <- lc_moments(df=df,orders=orders)
	}
	return(mu)
}
# cumulants of the log of the non-central chi-square
lnc_cumuls <- function(df,ncp=0,order.max=3,orders=c(1:order.max)) {
	mu <- lnc_moments(df,ncp,orders=orders)
	kappa <- PDQutils::moment2cumulant(mu)
	return(kappa)
}
# now proceed as usual! log noncentral chisquare!


# compute the cumulants of the sumlogchisq
# distribution. 
sumlogchisq_cumuls <- function(wts,df,ncp=0,order.max=3) {
	#nterms <- max(vapply(list(wts,df,ncp),length,0))
	subkappa <- mapply(function(w,dd,nn) 
										 { (w ^ (1:order.max)) * lnc_cumuls(df=dd,ncp=nn,order.max=order.max) },
										 wts,df,ncp,SIMPLIFY=FALSE)
	kappa <- Reduce('+', subkappa)
	return(kappa)
}

# dsumlogchisq, psumlogchisq, qsumlogchisq, rsumlogchisq#FOLDUP
#' @title The sum of the logs of (non-central) chi-squares distribution.
#'
#' @description 
#'
#' Density, distribution function, quantile function and random
#' generation for the distribution of the weighted sum of logs of
#' non-central chi-squares.
#'
#' @details
#'
#' Let \eqn{X_i \sim \chi^2\left(\delta_i, \nu_i\right)}{X_i ~ chi^2(delta_i, v_i)}
#' be independently distributed non-central chi-squares, where \eqn{\nu_i}{v_i}
#' are the degrees of freedom, and \eqn{\delta_i}{delta_i} are the
#' non-centrality parameters.  
#' Let \eqn{w_i} be given constants. Suppose
#' \deqn{Y = \sum_i w_i \log X_i.}{Y = sum w_i log(X_i).}
#' Then \eqn{Y}{Y} follows a weighted sum of log of chi-squares distribution. 
#'
#' @usage
#'
#' dsumlogchisq(x, wts, df, ncp=0, log = FALSE, order.max=6)
#'
#' psumlogchisq(q, wts, df, ncp=0, lower.tail = TRUE, log.p = FALSE, order.max=6)
#'
#' qsumlogchisq(p, wts, df, ncp=0, lower.tail = TRUE, log.p = FALSE, order.max=6)
#'
#' rsumlogchisq(n, wts, df, ncp=0)
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param wts the vector of weights. 
#' This is recycled against the \code{df, ncp}, but not against the \code{x,q,p,n}.
#' @param df the vector of degrees of freedom. 
#' This is recycled against the \code{wts, ncp}, but not against the \code{x,q,p,n}.
#' @param ncp the vector of non-centrality parameters. 
#' This is recycled against the \code{wts, df}, but not against the \code{x,q,p,n}.
#'
#' @template etc
#' @template distribution
#' @template apx_distribution
#' @template not-recycled
#' @template ref-lnc
#'
#' @return \code{dsumlogchisq} gives the density, \code{psumlogchisq} gives the 
#' distribution function, \code{qsumlogchisq} gives the quantile function, 
#' and \code{rsumlogchisq} generates random deviates.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @aliases dsumlogchisq psumlogchisq qsumlogchisq rsumlogchisq
#' @seealso 
#' The product of chi-squares to a power,
#' \code{\link{dprodchisqpow}},
#' \code{\link{pprodchisqpow}},
#' \code{\link{qprodchisqpow}},
#' \code{\link{rprodchisqpow}}.
#'
#' @examples 
#' wts <- c(1,-3,4)
#' df <- c(100,20,10)
#' ncp <- c(5,3,1)
#' rvs <- rsumlogchisq(128, wts, df, ncp)
#' dvs <- dsumlogchisq(rvs, wts, df, ncp)
#' qvs <- psumlogchisq(rvs, wts, df, ncp)
#' pvs <- qsumlogchisq(ppoints(length(rvs)), wts, df, ncp)
#' @rdname dsumlogchisq
#' @name sumlogchisq
#' @export 
dsumlogchisq <- function(x, wts, df, ncp=0, log = FALSE, order.max=6) {
	kappa <- sumlogchisq_cumuls(wts,df,ncp,order.max=order.max)
	retval <- PDQutils::dapx_edgeworth(x,kappa,log=log)
	return(retval)
}
#' @export
psumlogchisq <- function(q, wts, df, ncp=0, lower.tail = TRUE, log.p = FALSE, order.max=6) {
	kappa <- sumlogchisq_cumuls(wts,df,ncp,order.max=order.max)
	retval <- PDQutils::papx_edgeworth(q,kappa,lower.tail=lower.tail,log.p=log.p)
	return(retval)
}
#' @export 
qsumlogchisq <- function(p, wts, df, ncp=0, lower.tail = TRUE, log.p = FALSE, order.max=6) {
	kappa <- sumlogchisq_cumuls(wts,df,ncp,order.max=order.max)
	retval <- PDQutils::qapx_cf(p,kappa,lower.tail=lower.tail,log.p=log.p)
	return(retval)
}
#' @export 
rsumlogchisq <- function(n, wts, df, ncp=0) {
	subX <- mapply(function(w,dd,nn) { w * log(unbroken_rchisq(n,df=dd,ncp=nn)) },
										 wts,df,ncp,SIMPLIFY=FALSE)
	X <- Reduce('+', subX)
	return(X)
}
#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
