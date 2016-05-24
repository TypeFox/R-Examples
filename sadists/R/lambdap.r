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

# Created: 2014.02.21
# Copyright: Steven E. Pav, 2014
# Author: Steven E. Pav
# Comments: Steven E. Pav

# compute the cumulants of the lambda-prime
# distribution. this is distributed as
#
# t sqrt(chi^2(df) / df) + Z
#
# where the chi^2 is a chi-square independent of Z
lambdap_cumuls <- function(df,t,order.max=3) {
	# should check that they are scalars here
	stopifnot(length(df) == 1,length(t) == 1)
	kappa <- upsilon_cumuls(df,t,order.max)
	return(kappa)
}
# compute the raw moments of the lambda-prime
# distribution. this is distributed as
#
# t sqrt(chi^2(df) / df) + Z
#
# where the chi^2 is a chi-square independent of Z
lambdap_moms <- function(df,t,order.max=3) {
	mu <- cumulant2moment(lambdap_cumuls(df,t,order.max))
	return(mu)
}

# dlambdap, plambdap, qlambdap, rlambdap#FOLDUP
#' @title The lambda prime distribution.
#'
#' @description 
#'
#' Density, distribution function, quantile function and random
#' generation for the lambda prime distribution.
#'
#' @details
#'
#' Suppose \eqn{y \sim \chi^2\left(\nu\right)}{y ~ x^2(v)}, and
#' \eqn{Z}{Z} is a standard normal. 
#' \deqn{T = Z + t \sqrt{y/\nu}}{T = Z + t sqrt(y/v)}
#' takes a lambda prime distribution with parameters 
#' \eqn{\nu, t}{v, t}.
#' A lambda prime random variable can be viewed as a confidence
#' level on a non-central t because 
#' \deqn{t = \frac{Z' + T}{\sqrt{y/\nu}}}{t = (Z' + T)/sqrt(y/v)}
#'
#' @usage
#'
#' dlambdap(x, df, t, log = FALSE, order.max=6)
#'
#' plambdap(q, df, t, lower.tail = TRUE, log.p = FALSE, order.max=6)
#'
#' qlambdap(p, df, t, lower.tail = TRUE, log.p = FALSE, order.max=6)
#'
#' rlambdap(n, df, t)
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param df the degrees of freedom in the chi square. 
#' This is not recycled against the \code{x,q,p,n}.
#' @param t the scaling parameter on the chi.
#' This is not recycled against the \code{x,q,p,n}.
#'
#' @return \code{dlambdap} gives the density, \code{plambdap} gives the 
#' distribution function, \code{qlambdap} gives the quantile function, 
#' and \code{rlambdap} generates random deviates.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @aliases dlambdap plambdap qlambdap rlambdap
#' @seealso t distribution functions, \code{\link{dt}, \link{pt}, \link{qt}, \link{rt}},
#' K prime distribution functions, \code{\link{dkprime}, \link{pkprime}, \link{qkprime}, \link{rkprime}},
#' upsilon distribution functions, \code{\link{dupsilon}, \link{pupsilon}, \link{qupsilon}, \link{rupsilon}},
#' @template etc
#' @template ref-lambdap
#' @template distribution
#' @template apx_distribution
#' @template not-recycled
#' @examples 
#' rv <- rlambdap(100, 50, t=0.01)
#' d1 <- dlambdap(1, 50, t=0.01)
#' pv <- plambdap(rv, 50, t=0.01)
#' qv <- qlambdap(ppoints(length(rv)), 50, t=1)
#'
#' @name lambdap
#' @rdname dlambdap
#' @export 
dlambdap <- function(x, df, t, log = FALSE, order.max=6) {
	kappa <- lambdap_cumuls(df,t,order.max=order.max)
	retval <- PDQutils::dapx_edgeworth(x,kappa,log=log)
	return(retval)
}
#' @export 
plambdap <- function(q, df, t, lower.tail = TRUE, log.p = FALSE, order.max=6) {
	kappa <- lambdap_cumuls(df,t,order.max=order.max)
	retval <- PDQutils::papx_edgeworth(q,kappa,lower.tail=lower.tail,log.p=log.p)
	return(retval)
}
#' @export 
qlambdap <- function(p, df, t, lower.tail = TRUE, log.p = FALSE, order.max=6) {
	kappa <- lambdap_cumuls(df,t,order.max=order.max)
	retval <- PDQutils::qapx_cf(p,kappa,lower.tail=lower.tail,log.p=log.p)
	return(retval)
}
#' @export 
rlambdap <- function(n, df, t) {
	y <- rchisq(n,df=df) 
	X <- rnorm(n,mean=t * sqrt(y/df))
	return(X)
}
#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
