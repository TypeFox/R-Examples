# Copyright 2012-2014 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav

# This file is part of SharpeR.
#
# SharpeR is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SharpeR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with SharpeR.  If not, see <http://www.gnu.org/licenses/>.

# env var:
# nb: 
# see also:
# todo:
# changelog: 
#
# Created: 2012.05.19
# Copyright: Steven E. Pav, 2012-2013
# Author: Steven E. Pav
# Comments: Steven E. Pav

#' @include utils.r

# 2FIX: is df the # of observations or the d.f. of the t-stat? ack!

# note: on citations, use the Chicago style from google scholar. tks.

# the following are useful for grokking R
# showMethods("print") 
# getAnywhere("t.test.default")

# roxygen2 tiparoonie on S3 methods:
# http://stackoverflow.com/a/7199577/164611

# @param log,log.p logical; if TRUE, probabilities p are given as \eqn{\mbox{log}(p)}{log(p)}.
# @param lower.tail logical; if TRUE (default), probabilities are
#        \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}{P[X > x]}.

########################################################################
# Distributions
########################################################################

# some facts about t-stats#FOLDUP
# geometric bias in the expectation of the non-central t-stat with a
# given number of d.f.
# this is 1 over 'c4', essentially:
# http://mathworld.wolfram.com/StandardDeviationDistribution.html
# http://finzi.psych.upenn.edu/R/library/IQCC/html/c4.html
# http://math.stackexchange.com/questions/71573/the-limit-of-the-ratio-of-two-gammax-functions
#.tbias <- function(df) { 
	#retv <- sqrt(df / 2) * exp(lgamma((df-1)/2) - lgamma(df/2))
	#return(retv)
#}
# same, but for sr:
#.srbias <- function(df) { 
	#return(.tbias(df-1))
#}
#UNFOLD

# rescaled t-distributions
# drt, prt, qrt, rrt#FOLDUP
#
# the rescaled t-distributions. r follows a rescaled t-distribution
# with m df, rescaling K, and non-centrality parameter rho if
# f/K follows a non-central t with m df and non-centrality parameter rho.
# this is a very thin wrapper.
drt <- function(x, df, K, rho = 0, log = FALSE) {
	tx <- x / K
	if (missing(rho)) {
		xd <- dt(tx, df = df, log = log)
	} else {
		ncp <- rho / K
		xd <- dt(tx, df = df, ncp = ncp, log = log)
	}
	if (log) {
		retv <- xd - log(K)
	} else {
		retv <- xd / K
	}
	return(retv)	
}
prt <- function(q, df, K, rho = 0, ...) {
	tq <- q / K
	if (missing(rho)) {
		retv <- pt(q = tq, df = df, ...)		
	} else {
		ncp <- rho / K
		retv <- pt(q = tq, df = df, ncp = ncp, ...)		
	}
	return(retv)	
}
qrt <- function(p, df, K, rho = 0, ...) {
	if (missing(rho)) {
		tq <- qt(p, df = df, ...)
	} else {
		ncp <- rho / K
		tq <- qt(p, df = df, ncp, ...)
	}
	retv <- tq * K
	return(retv)
}
rrt <- function(n, df, K, rho = 0) {
	if (missing(rho)) {
		tr <- rt(n, df = df)
	} else {
		ncp <- rho / K
		tr <- rt(n, df = df, ncp = ncp)
	}
	retv <- tr * K
	return(retv)
}
#UNFOLD

# Sharpe ratio as a distribution
# dsr, psr, qsr, rsr#FOLDUP
#' @title The (non-central) Sharpe ratio.
#'
#' @description 
#'
#' Density, distribution function, quantile function and random
#' generation for the Sharpe ratio distribution with \code{df} degrees of freedom
#' (and optional signal-noise-ratio \code{zeta}).
#'
#' @details
#'
#' Suppose \eqn{x_i}{xi} are \eqn{n} independent draws of a normal random
#' variable with mean \eqn{\mu}{mu} and variance \eqn{\sigma^2}{sigma^2}.
#' Let \eqn{\bar{x}}{xbar} be the sample mean, and \eqn{s} be
#' the sample standard deviation (using Bessel's correction). Let \eqn{c_0}{c0}
#' be the 'risk free rate'.  Then
#' \deqn{z = \frac{\bar{x} - c_0}{s}}{z = (xbar - c0)/s} 
#' is the (sample) Sharpe ratio.
#' 
#' The units of \eqn{z} is \eqn{\mbox{time}^{-1/2}}{per root time}.
#' Typically the Sharpe ratio is \emph{annualized} by multiplying by
#' \eqn{\sqrt{d}}{sqrt(d)}, where \eqn{d} is the number of observations
#' per epoch (typically a year).
#'
#' Letting \eqn{z = \sqrt{d}\frac{\bar{x}-c_0}{s}}{z = sqrt(d)(xbar - c0)/s},
#' where the sample estimates are based on \eqn{n} observations, 
#' then \eqn{z}{z} takes a (non-central) Sharpe ratio distribution
#' parametrized by \eqn{n} 'degrees of freedom', non-centrality parameter
#' \eqn{\zeta = \frac{\mu - c_0}{\sigma}}{zeta = (mu - c0)/sigma}, and 
#' annualization parameter \eqn{d}. 
#'
#' The parameters are encoded as follows:
#' \itemize{
#' \item \eqn{n} is denoted by \code{df}.
#' \item \eqn{\zeta}{zeta} is denoted by \code{zeta}.
#' \item \eqn{d} is denoted by \code{ope}. ('Observations Per Year')
#' }
#' 
#' If the returns violate the assumptions of normality, independence, etc
#' (\emph{as they always should in the real world}), the sample Sharpe Ratio
#' will not follow this distribution. It does provide, however, a reasonable
#' approximation in many cases.
#'
#' @usage
#'
#' dsr(x, df, zeta, ope, ...)
#'
#' psr(q, df, zeta, ope, ...)
#'
#' qsr(p, df, zeta, ope, ...)
#'
#' rsr(n, df, zeta, ope)
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param df the number of observations the statistic is based on. This 
#'        is one more than the number of degrees of freedom in the
#'        corresponding t-statistic, although the effect will be small
#'        when \code{df} is large.
#' @param zeta the 'signal-to-noise' parameter, \eqn{\zeta}{zeta} defined as the population
#'        mean divided by the population standard deviation, 'annualized'.
#' @template param-ope
#' @param ... arguments passed on to the respective t-distribution functions, namely
#' \code{lower.tail} with default \code{TRUE}, \code{log} with default \code{FALSE}, 
#' and \code{log.p} with default \code{FALSE}.
#' @keywords distribution 
#' @return \code{dsr} gives the density, \code{psr} gives the distribution function,
#' \code{qsr} gives the quantile function, and \code{rsr} generates random deviates.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @rdname dsr
#' @aliases psr qsr rsr
#' @seealso t-distribution functions, \code{\link{dt}, \link{pt}, \link{qt}, \link{rt}}
#' @note
#' This is a thin wrapper on the t distribution. 
#' The functions \code{\link{dt}, \link{pt}, \link{qt}} can accept ncp from
#' limited range (\eqn{|\delta|\le 37.62}{delta <= 37.62}). Some corrections
#' may have to be made here for large \code{zeta}.
#' @export 
#' @template etc
#' @template sr
#' @examples 
#' rvs <- rsr(128, 253*6, 0, 253)
#' dvs <- dsr(rvs, 253*6, 0, 253)
#' pvs.H0 <- psr(rvs, 253*6, 0, 253)
#' pvs.HA <- psr(rvs, 253*6, 1, 253)
#' \dontrun{
#' plot(ecdf(pvs.H0))
#' plot(ecdf(pvs.HA))
#' }
#'
dsr <- function(x, df, zeta, ope = 1, ...) {
	K <- sqrt(ope / df)
	if (missing(zeta)) {
		retv <- drt(x, df-1, K, ...)
	} else {
		retv <- drt(x, df-1, K, rho = zeta, ...)
	}
	return(retv)	
}
#' @export 
psr <- function(q, df, zeta, ope, ...) {
	K <- sqrt(ope / df)
	if (missing(zeta)) {
		retv <- prt(q, df-1, K, ...)
	} else {
		retv <- prt(q, df-1, K, rho = zeta, ...)
	}
	return(retv)	
}
#' @export 
qsr <- function(p, df, zeta, ope, ...) {
	K <- sqrt(ope / df)
	if (missing(zeta)) {
		retv <- qrt(p, df-1, K, ...)
	} else {
		retv <- qrt(p, df-1, K, rho = zeta, ...)
	}
	return(retv)
}
#' @export 
rsr <- function(n, df, zeta, ope) {
	K <- sqrt(ope / df)
	if (missing(zeta)) {
		retv <- rrt(n, df-1, K)
	} else {
		retv <- rrt(n, df-1, K, rho = zeta)
	}
	return(retv)
}
#UNFOLD

# Hotelling
# dT2, pT2, qT2, rT2#FOLDUP
#  ' @title The (non-central) Hotelling distribution.
#  '
#  ' @description 
#  '
#  ' Density, distribution function, quantile function and random
#  ' generation for the Hotelling distribution distribution with 
#  ' \code{df1} and \code{df2} degrees of freedom
#  ' (and optional non-centrality parameter \code{delta2}).
#  '
#  ' @details
#  '
#  ' Suppose \eqn{x_i}{xi} are \eqn{n}{n} independent draws of a \eqn{q}{q}-variate
#  ' normal random variable with mean \eqn{\mu}{mu} and covariance matrix
#  ' \eqn{\Sigma}{Sigma}. Let \eqn{\bar{x}}{xbar} be the (vector) sample mean, and 
#  ' \eqn{S}{S} be the sample covariance matrix (using Bessel's correction). Then
#  ' \deqn{T^2 = n \bar{x}^{\top}S^{-1}\bar{x}}{T^2 = n xbar' S^-1 xbar} 
#  ' follows a (non-central)
#  ' Hotelling T-squared distribution with \eqn{q}{q} and \eqn{n-1}{n-1}
#  ' degrees of freedom, and non-centrality parameter
#  ' \deqn{\delta^2 = n \mu^{\top}\Sigma^{-1}\mu}{delta^2 = n mu' Sigma^-1 mu}
#  '
#  ' The (non-central) T-squared distribution is a (non-central) F distribution
#  ' up to scaling which depends on \eqn{q}{q} and \eqn{n}{n}.
#  '
#  ' @usage
#  '
#  ' dT2(x, df1, df2, delta2, log = FALSE)
#  '
#  ' pT2(q, df1, df2, delta2, ...)
#  '
#  ' qT2(p, df1, df2, delta2, ...)
#  '
#  ' rT2(n, df1, df2, delta2)
#  '
#  ' @param x,q vector of quantiles.
#  ' @param p vector of probabilities.
#  ' @param n number of observations. 
#  ' @param df1 the dimension of the vector space from which multivariate
#  '        observations had been drawn, \eqn{q}{q}.
#  ' @param df2 the number of observations that the sample mean and covariance
#  '        are based on, \eqn{n}{n}.
#  ' @param delta2 the population non-centrality parameter, defined as 
#  '        \eqn{\delta^2 = n \mu^{\top}\Sigma^{-1}\mu}{delta^2 = n (mu' Sigma^-1 mu)}.
#  '        defaults to 0, i.e. a central \eqn{T^2}{T2} distribution.
#  ' @param log logical; if TRUE, probabilities p are given as \eqn{\mbox{log}(p)}{log(p)}.
#  ' @param ... arguments passed on to the respective F distribution functions, namely
#  ' \code{lower.tail} with default \code{TRUE}, and \code{log.p}, with default \code{FALSE}.
#  ' @keywords distribution 
#  ' @return \code{dT2} gives the density, \code{pT2} gives the distribution function,
#  ' \code{qT2} gives the quantile function, and \code{rT2} generates random deviates.
#  '
#  ' Invalid arguments will result in return value \code{NaN} with a warning.
#  ' @aliases pT2 
#  ' @aliases qT2 
#  ' @aliases rT2 
#  ' @seealso F-distribution functions, \code{\link{df}, \link{pf}, \link{qf}, \link{rf}},
#  ' Sharpe ratio distribution, \code{\link{dsr}, \link{psr}, \link{qsr}, \link{rsr}}.
#  ' @template etc
#  ' @family Hotelling
#  ' @note
#  ' This is a thin wrapper on the F distribution, provided for convenience.
#  ' @references 
#  '
#  ' Wikipedia contributors, 'Hotelling's T-squared distribution', Wikipedia, The Free Encyclopedia, 
#  ' 11 December 2012, 13:38 UTC, \url{http://en.wikipedia.org/w/index.php?title=Hotelling\%27s_T-squared_distribution\&oldid=527530524}
#  ' [accessed 9 January 2013]
#  '
#  ' Bilodeau, Martin, and David Brenner. Theory of multivariate statistics. Springer, 1999.
#  ' \url{http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.172.3290}
#  '
#  ' Timm, Neil H. Applied multivariate analysis: methods and case studies. Springer, 2002.
#  ' \url{http://books.google.com/books?id=vtiyg6fnnskC}
#  ' 
#  ' Hotelling, Harold. "The Generalization of Student's Ratio." Annals of Mathematical 
#  ' Statistics 2, no. 3 (1931): 360--378. \url{http://projecteuclid.org/euclid.aoms/1177732979}
#  '
#  ' @examples 
#  ' rvs <- rT2(128, 4, 253*6, 0)
#  ' dvs <- dT2(rvs, 4, 253*6, 0)
#  ' pvs <- pT2(rvs, 4, 253*6, 0)
#  ' plot(ecdf(pvs))
#  ' pvs <- pT2(rvs, 4, 253*6, 1)
#  ' plot(ecdf(pvs))
#  '
#  ' @export
dT2 <- function(x, df1, df2, delta2, log = FALSE) {
	Fs <- .T2_to_F(x, df1, df2)
	if (missing(delta2)) {
		dv <- df(Fs, df1 = df1, df2 = df2 - df1, log = log)
	} else {
		dv <- df(Fs, df1 = df1, df2 = df2 - df1, ncp = delta2, log = log)
	}
	if (log) {
		retv <- (dv + log(.d_T2_to_F(x, df1, df2)))
	} else {
		retv <- (dv * .d_T2_to_F(x, df1, df2))
	}
	return(retv)
}
#  ' @export
pT2 <- function(q, df1, df2, delta2, ...) {
	Fs <- .T2_to_F(q, df1, df2)
	#cat(sprintf('pT2 F(%d,%d) stat: %g\n',df1,df2,Fs))
	if (missing(delta2)) {
		#retv <- pf(Fs, df1 = df1, df2 = df2 - df1, ncp = 0, ...)
		retv <- pf(Fs, df1 = df1, df2 = df2 - df1, ...)
	} else {
		retv <- pf(Fs, df1 = df1, df2 = df2 - df1, ncp = delta2, ...)
	}
	return(retv)
}
#  ' @export
qT2 <- function(p, df1, df2, delta2, ... ) {
	if (missing(delta2)) {
		#Fq <- qf(Fp, df1 = df1, df2 = df2 - df1, ncp = 0, ... )
		Fq <- qf(p, df1 = df1, df2 = df2 - df1, ... )
	} else {
		Fq <- qf(p, df1 = df1, df2 = df2 - df1, ncp = delta2, ... )
	}
	retv <- .F_to_T2(Fq, df1, df2)
	return(retv)
}
#  ' @export
rT2 <- function(n, df1, df2, delta2) {
	if (missing(delta2)) {
		Fr <- rf(n, df1 = df1, df2 = df2 - df1)
	} else {
		Fr <- rf(n, df1 = df1, df2 = df2 - df1, ncp = delta2)
	}
	retv <- .F_to_T2(Fr, df1, df2)
	return(retv)
}
#UNFOLD

# SR^*
# dsropt, psropt, qsropt, rsropt#FOLDUP
#' @title The (non-central) maximal Sharpe ratio distribution.
#'
#' @description 
#'
#' Density, distribution function, quantile function and random
#' generation for the maximal Sharpe ratio distribution with 
#' \code{df1} and \code{df2} degrees of freedom
#' (and optional maximal signal-noise-ratio \code{zeta.s}).
#'
#' @details
#'
#' Suppose \eqn{x_i}{xi} are \eqn{n} independent draws of a \eqn{q}-variate
#' normal random variable with mean \eqn{\mu}{mu} and covariance matrix
#' \eqn{\Sigma}{Sigma}. Let \eqn{\bar{x}}{xbar} be the (vector) sample mean, and 
#' \eqn{S} be the sample covariance matrix (using Bessel's correction). Let
#' \deqn{Z(w) = \frac{w^{\top}\bar{x} - c_0}{\sqrt{w^{\top}S w}}}{Z(w) = (w'xbar - c0)/sqrt(w'Sw)}
#' be the (sample) Sharpe ratio of the portfolio \eqn{w}, subject to 
#' risk free rate \eqn{c_0}{c0}.
#'
#' Let \eqn{w_*}{w*} be the solution to the portfolio optimization problem:
#' \deqn{\max_{w: 0 < w^{\top}S w \le R^2} Z(w),}{max {Z(w) | 0 < w'Sw <= R^2},}
#' with maximum value \eqn{z_* = Z\left(w_*\right)}{z* = Z(w*)}.
#' Then 
#' \deqn{w_* = R \frac{S^{-1}\bar{x}}{\sqrt{\bar{x}^{\top}S^{-1}\bar{x}}}}{%
#' w* = R S^-1 xbar / sqrt(xbar' S^-1 xbar)}
#' and
#' \deqn{z_* = \sqrt{\bar{x}^{\top} S^{-1} \bar{x}} - \frac{c_0}{R}}{%
#' z* = sqrt(xbar' S^-1 xbar) - c0/R}
#'
#' The variable \eqn{z_*}{z*} follows an \emph{Optimal Sharpe ratio}
#' distribution. For convenience, we may assume that the sample statistic
#' has been annualized in the same manner as the Sharpe ratio, that is 
#' by multiplying by \eqn{d}, the number of observations per
#' epoch.
#' 
#' The Optimal Sharpe Ratio distribution is parametrized by the number 
#' of assets, \eqn{q}, the number of independent observations, \eqn{n}, the 
#' noncentrality parameter, 
#' \deqn{\zeta_* = \sqrt{\mu^{\top}\Sigma^{-1}\mu},}{zeta* = sqrt(mu' Sigma^-1 mu),}
#' the 'drag' term, \eqn{c_0/R}{c0/R}, and the annualization factor, \eqn{d}.
#' The drag term makes this a location family of distributions, and 
#' by default we assume it is zero.
#' 
#' The parameters are encoded as follows:
#' \itemize{
#' \item \eqn{q} is denoted by \code{df1}.
#' \item \eqn{n} is denoted by \code{df2}.
#' \item \eqn{\zeta_*}{zeta*} is denoted by \code{zeta.s}.
#' \item \eqn{d} is denoted by \code{ope}.
#' \item \eqn{c_0/R} is denoted by \code{drag}.
#' }
#'
#' @usage
#'
#' dsropt(x, df1, df2, zeta.s, ope, drag = 0, log = FALSE)
#'
#' psropt(q, df1, df2, zeta.s, ope, drag = 0, ...)
#'
#' qsropt(p, df1, df2, zeta.s, ope, drag = 0, ...)
#'
#' rsropt(n, df1, df2, zeta.s, ope, drag = 0, ...)
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. 
#' @param df1 the number of assets in the portfolio.
#' @param df2 the number of observations.
#' @param zeta.s the non-centrality parameter, defined as 
#'        \eqn{\zeta_* = \sqrt{\mu^{\top}\Sigma^{-1}\mu},}{zeta* = sqrt(mu' Sigma^-1 mu),}
#'        for population parameters.
#'        defaults to 0, \emph{i.e.} a central maximal Sharpe ratio distribution.
#' @template param-ope
#' @param drag the 'drag' term, \eqn{c_0/R}{c0/R}. defaults to 0. It is assumed
#'        that \code{drag} has been annualized, \emph{i.e.} is given in the
#'        same units as \code{x} and \code{q}.
#' @param log logical; if TRUE, densities \eqn{f} are given as \eqn{\mbox{log}(f)}{log(f)}.
#' @param ... arguments passed on to the respective Hotelling \eqn{T^2} functions.
#' @keywords distribution 
#' @return \code{dsropt} gives the density, \code{psropt} gives the distribution function,
#' \code{qsropt} gives the quantile function, and \code{rsropt} generates random deviates.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @rdname dsropt
#' @aliases psropt qsropt rsropt
#' @seealso F-distribution functions, \code{\link{df}, \link{pf}, \link{qf}, \link{rf}}, 
#' Sharpe ratio distribution, \code{\link{dsr}, \link{psr}, \link{qsr}, \link{rsr}}.
#' @export 
#' @template etc
#' @template sropt
#' @note
#' This is a thin wrapper on the Hotelling T-squared distribution, which is a
#' wrapper on the F distribution.
#' @references 
#'
#' Kan, Raymond and Smith, Daniel R. "The Distribution of the Sample Minimum-Variance Frontier."
#' Journal of Management Science 54, no. 7 (2008): 1364--1380.
#' \url{http://mansci.journal.informs.org/cgi/doi/10.1287/mnsc.1070.0852}
#'
#' @examples 
#' # generate some variates 
#' ngen <- 128
#' ope <- 253
#' df1 <- 8
#' df2 <- ope * 10
#' drag <- 0
#' # sample
#' rvs <- rsropt(ngen, df1, df2, drag, ope)
#' hist(rvs)
#' # these should be uniform:
#' isp <- psropt(rvs, df1, df2, drag, ope)
#' plot(ecdf(isp))
#'
dsropt <- function(x, df1, df2, zeta.s, ope, drag = 0, log = FALSE) {
	if (!missing(drag) && (drag != 0)) {
		x <- x + drag
	}
	if (!missing(ope)) {
		x <- .deannualize(x, ope)
		if (!missing(zeta.s)) {
			zeta.s <- .deannualize(zeta.s, ope)
		}
	}
	x.T2 <- .sropt_to_T2(x, df2)
	if (missing(zeta.s)) {
		delta2 <- 0
	} else {
		delta2 <- .sropt_to_T2(zeta.s, df2)
	}
	d.T2 <- dT2(x.T2, df1, df2, delta2, log=log)
	if (log) {
		retv <- (d.T2 - log(.d_T2_to_sropt(x, df2)))
	} else {
		retv <- (d.T2 / .d_T2_to_sropt(x, df2))
	}
	return(retv)
}
#' @export 
psropt <- function(q, df1, df2, zeta.s, ope, drag = 0, ...) {
	if (!missing(drag) && (drag != 0)) {
		q <- q + drag
	}
	if (!missing(ope)) {
		q <- .deannualize(q, ope)
		if (!missing(zeta.s)) {
			zeta.s <- .deannualize(zeta.s, ope)
		}
	}
	q.T2 <- .sropt_to_T2(q, df2)
	if (missing(zeta.s)) {
		delta2 = 0.0
	} else {
		delta2 <- .sropt_to_T2(zeta.s, df2)
	}
	retv <- pT2(q.T2, df1, df2, delta2, ...)
	return(retv)
}
#' @export 
qsropt <- function(p, df1, df2, zeta.s, ope, drag = 0, ...) {
	if (missing(zeta.s)) {
		delta2 = 0.0
	} else {
		if (!missing(ope)) {
			zeta.s <- .deannualize(zeta.s, ope)
		}
		delta2 <- .sropt_to_T2(zeta.s, df2)
	}
	q.T2 <- qT2(p, df1, df2, delta2, ...)
	retv <- .T2_to_sropt(q.T2, df2)
	if (!missing(ope)) {
		retv <- .annualize(retv,ope)
	}
	if (!missing(drag) && (drag != 0)) {
		retv <- retv - drag
	}
	return(retv)
}
#' @export 
rsropt <- function(n, df1, df2, zeta.s, ope, drag = 0, ...) {
	if (missing(zeta.s)) {
		delta2 = 0.0
	} else {
		if (!missing(ope)) {
			zeta.s <- .deannualize(zeta.s, ope)
		}
		delta2 <- .sropt_to_T2(zeta.s, df2)
	}
	r.T2 <- rT2(n, df1, df2, delta2, ...) 
	retv <- .T2_to_sropt(r.T2, df2)
	if (!missing(ope)) {
		retv <- .annualize(retv,ope)
	}
	if (!missing(drag) && (drag != 0)) {
		retv <- retv - drag
	}
	return(retv)
}
#UNFOLD

########################################################################
# 'confidence distributions'

# lambda prime
# plambdap, qlambdap, rlambdap#FOLDUP
#' @title The lambda-prime distribution.
#'
#' @description 
#'
#' Distribution function and quantile function for LeCoutre's
#' lambda-prime distribution with \code{df} degrees of freedom
#' and the observed t-statistic, \code{tstat}.
#'
#' @details
#'
#' Let \eqn{t}{t} be distributed
#' as a non-central t with \eqn{\nu}{v} degrees of freedom and non-centrality
#' parameter \eqn{\delta}{ncp}. We can view this as
#' \deqn{t = \frac{Z + \delta}{\sqrt{V/\nu}}.}{t = (Z + ncp)/sqrt(V/v)}
#' where \eqn{Z}{Z} is a standard normal, \eqn{\delta}{ncp} is the
#' non-centrality parameter, \eqn{V}{V} is a chi-square RV with \eqn{\nu}{v}
#' degrees of freedom, independent of \eqn{Z}{Z}.  We can rewrite this as
#' \deqn{\delta = t\sqrt{V/\nu} + Z.}{ncp = t sqrt(V/v) + Z}
#' 
#' Thus a 'lambda-prime' random variable with parameters \eqn{t}{t} and
#' \eqn{\nu}{v} is one expressable as a sum
#' \deqn{t\sqrt{V/\nu} + Z}{t sqrt(V/v) + Z}
#' for Chi-square \eqn{V}{V} with \eqn{\nu}{v} d.f., independent from
#' standard normal \eqn{Z}{Z}
#'
#' @usage
#'
#' plambdap(q, df, tstat, lower.tail = TRUE, log.p = FALSE)
#'
#' qlambdap(p, df, tstat, lower.tail = TRUE, log.p = FALSE)
#'
#' rlambdap(n, df, tstat)
#'
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If 'length(n) > 1', the length is
#' taken to be the number required.
#' @param df the degrees of freedom of the t-statistic.
#' @param tstat the observed (non-central) t-statistic.
#' @param log.p logical; if TRUE, probabilities p are given as \eqn{\mbox{log}(p)}{log(p)}.
#' @param lower.tail logical; if TRUE (default), probabilities are
#'        \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}{P[X > x]}.
#' @inheritParams dsr
#' @keywords distribution 
#' @return \code{dlambdap} gives the density, \code{plambdap} gives the distribution function,
#' \code{qlambdap} gives the quantile function, and \code{rlambdap} generates random deviates.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @aliases plambdap qlambdap rlambdap
#' @seealso t-distribution functions, \code{\link{dt},\link{pt},\link{qt},\link{rt}}
#' @export 
#' @template etc
#' @family sr
#' @references 
#'
#' Lecoutre, Bruno. "Another look at confidence intervals for the noncentral t distribution." 
#' Journal of Modern Applied Statistical Methods 6, no. 1 (2007): 107--116.
#' \url{http://www.univ-rouen.fr/LMRS/Persopage/Lecoutre/telechargements/Lecoutre_Another_look-JMSAM2007_6(1).pdf}
#'
#' Lecoutre, Bruno. "Two useful distributions for Bayesian predictive procedures under normal models."
#' Journal of Statistical Planning and Inference 79  (1999): 93--105. 
#'
#' @note
#' \code{plambdap} should be an increasing function of the argument \code{q},
#' and decreasing in \code{tstat}. \code{qlambdap} should be increasing
#' in \code{p}
#' @examples 
#' rvs <- rnorm(128)
#' pvs <- plambdap(rvs, 253*6, 0.5)
#' plot(ecdf(pvs))
#' pvs <- plambdap(rvs, 253*6, 1)
#' plot(ecdf(pvs))
#' pvs <- plambdap(rvs, 253*6, -0.5)
#' plot(ecdf(pvs))
#' # test vectorization:
#' qv <- qlambdap(0.1,128,2)
#' qv <- qlambdap(c(0.1),128,2)
#' qv <- qlambdap(c(0.2),128,2)
#' qv <- qlambdap(c(0.2),253,2)
#' qv <- qlambdap(c(0.1,0.2),128,2)
#' qv <- qlambdap(c(0.1,0.2),c(128,253),2)
#' qv <- qlambdap(c(0.1,0.2),c(128,253),c(2,4))
#' qv <- qlambdap(c(0.1,0.2),c(128,253),c(2,4,8,16))
#' # random generation
#' rv <- rlambdap(1000,252,2)
#'
plambdap <- function(q,df,tstat,lower.tail=TRUE,log.p=FALSE) {
	# this is just a silly wrapper on pt
	retv <- pt(q=tstat,df=df,ncp=q,lower.tail=!lower.tail,log.p=log.p)
	return(retv)
}

# create a scalar function that we later vectorize. 
.qlambdap <- function(p,df,tstat,lower.tail=TRUE,log.p=FALSE) { # nocov start
	if (!log.p) {
		if (p == 1)
			return(ifelse(lower.tail,Inf,-Inf))
		if (p == 0)
			return(ifelse(lower.tail,-Inf,Inf))
		if ((p < 0) || (p > 1))
			return (NaN)
	} else {
		if (p == 0)
			return(ifelse(lower.tail,Inf,-Inf))
		if (p > 0)
			return (NaN)
		if (is.infinite(p))
			return(ifelse(lower.tail,-Inf,Inf))
	}

	# create a function increasing in its argument that
	# we wish to zero
	if (lower.tail) {
		zerf <- function(q) {
			return(plambdap(q,df=df,tstat=tstat,lower.tail=lower.tail,log.p=log.p) - p)
		}
	} else {
		zerf <- function(q) {
			return(p - plambdap(q,df=df,tstat=tstat,lower.tail=lower.tail,log.p=log.p))
		}
	}
	zmax = 2 * max(qnorm(p,lower.tail=TRUE,log.p=log.p),qnorm(p,lower.tail=FALSE,log.p=log.p))
	flo <- min(-1,tstat - zmax)
	fhi <- max( 1,tstat + zmax)

	#find approximate endpoints;
	#expand them until pt(tstat,df,flo) > p and pt(tstat,df,fhi) < p
	FLIM <- 1e8
	while ((zerf(flo) > 0) && (flo > -FLIM)) { flo <- 2 * flo }
	while ((zerf(fhi) < 0) && (fhi < FLIM))  { fhi <- 2 * fhi }

	ncp <- uniroot(zerf,c(flo,fhi))
	return(ncp$root)
} # nocov end
#' @export 
qlambdap <- Vectorize(.qlambdap, 
											vectorize.args = c("p","df","tstat"),
											SIMPLIFY = TRUE)
#' @export 
rlambdap <- function(n, df, tstat) {
	rvs <- rnorm(n) + tstat * sqrt(rchisq(n, df=df) / df)
}
#UNFOLD

# co-SR^*
# pco_sropt, qco_sropt#FOLDUP
#' @title The 'confidence distribution' for maximal Sharpe ratio.
#'
#' @description 
#'
#' Distribution function and quantile function for the 'confidence
#' distribution' of the maximal Sharpe ratio. This is just an inversion
#' to perform inference on \eqn{\zeta_*}{zeta*} given observed statistic 
#' \eqn{z_*}{z*}.
#'
#' @details
#' 
#' Suppose \eqn{z_*}{z*} follows a \emph{Maximal Sharpe ratio} distribution
#' (see \code{\link{SharpeR-package}}) for known degrees of freedom, and 
#' unknown non-centrality parameter \eqn{\zeta_*}{zeta*}. The 
#' 'confidence distribution' views \eqn{\zeta_*}{zeta*} as a random
#' quantity once \eqn{z_*}{z*} is observed. As such, the CDF of
#' the confidence distribution is the same as that of the 
#' Maximal Sharpe ratio (up to a flip of \code{lower.tail});
#' while the quantile function is used to compute confidence
#' intervals on \eqn{\zeta_*}{zeta*} given \eqn{z_*}{z*}.
#'
#' @usage
#'
#' pco_sropt(q,df1,df2,z.s,ope,lower.tail=TRUE,log.p=FALSE) 
#'
#' qco_sropt(p,df1,df2,z.s,ope,lower.tail=TRUE,log.p=FALSE,lb=0,ub=Inf) 
#'
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param z.s an observed Sharpe ratio statistic, annualized.
#' @template param-ope
#' @param log.p logical; if TRUE, probabilities p are given as \eqn{\mbox{log}(p)}{log(p)}.
#' @param lower.tail logical; if TRUE (default), probabilities are
#'        \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}{P[X > x]}.
#' @param lb the lower bound for the output of \code{qco_sropt}.
#' @param ub the upper bound for the output of \code{qco_sropt}.
#' @inheritParams dsropt
#' @inheritParams dsr
#' @keywords distribution 
#' @return \code{pco_sropt} gives the distribution function, and
#' \code{qco_sropt} gives the quantile function.
#'
#' Invalid arguments will result in return value \code{NaN} with a warning.
#' @aliases qco_sropt 
#' @seealso \code{\link{dsropt},\link{psropt},\link{qsropt},\link{rsropt}}
#' @export 
#' @template etc
#' @family sropt
#' @note
#' When \code{lower.tail} is true, \code{pco_sropt} is monotonic increasing 
#' with respect to \code{q}, and decreasing in \code{sropt}; these are reversed
#' when \code{lower.tail} is false. Similarly, \code{qco_sropt} is increasing
#' in \code{sign(as.double(lower.tail) - 0.5) * p} and
#' \code{- sign(as.double(lower.tail) - 0.5) * sropt}.
#'
#' @examples 
#'
#' zeta.s <- 2.0
#' ope <- 253
#' ntest <- 50
#' df1 <- 4
#' df2 <- 6 * ope
#' rvs <- rsropt(ntest,df1=df1,df2=df2,zeta.s=zeta.s)
#' qvs <- seq(0,10,length.out=51)
#' pps <- pco_sropt(qvs,df1,df2,rvs[1],ope)
#' \dontrun{
#' if (require(txtplot))
#'  txtplot(qvs,pps)
#' }
#' pps <- pco_sropt(qvs,df1,df2,rvs[1],ope,lower.tail=FALSE)
#' \dontrun{
#' if (require(txtplot))
#'  txtplot(qvs,pps)
#' }
#' 
#' svs <- seq(0,4,length.out=51)
#' pps <- pco_sropt(2,df1,df2,svs,ope)
#' pps <- pco_sropt(2,df1,df2,svs,ope,lower.tail=FALSE)
#' 
#' pps <- pco_sropt(qvs,df1,df2,rvs[1],ope,lower.tail=FALSE)
#' pco_sropt(-1,df1,df2,rvs[1],ope)
#'
#' qvs <- qco_sropt(0.05,df1=df1,df2=df2,z.s=rvs)
#' mean(qvs > zeta.s)
#' qvs <- qco_sropt(0.5,df1=df1,df2=df2,z.s=rvs)
#' mean(qvs > zeta.s)
#' qvs <- qco_sropt(0.95,df1=df1,df2=df2,z.s=rvs)
#' mean(qvs > zeta.s)
#' # test vectorization:
#' qv <- qco_sropt(0.1,df1,df2,rvs)
#' qv <- qco_sropt(c(0.1,0.2),df1,df2,rvs)
#' qv <- qco_sropt(c(0.1,0.2),c(df1,2*df1),df2,rvs)
#' qv <- qco_sropt(c(0.1,0.2),c(df1,2*df1),c(df2,2*df2),rvs)
#'
pco_sropt <- function(q,df1,df2,z.s,ope=1,lower.tail=TRUE,log.p=FALSE) {
	# 2FIX: add ope?
	# 2FIX: do the annualization just once for efficiency?
	# this is just a silly wrapper on psropt
	# delegate
	retv <- psropt(q=z.s,df1=df1,df2=df2,zeta.s=q,ope=ope,
									lower.tail=!lower.tail,log.p=log.p)  # sic the tail reversal
	return(retv)
}
# create a scalar function that we later vectorize. 
# 
# this inverts pco_sropt; note that when lower.tail=TRUE,
# pco_sropt is increasing in q, but decreasing in sropt.
# pco_sropt only accepts non-negative q 
#
# here we try to find lb <= q < ub such that
# pco_sropt(q,df1,df2,sropt,ope,lower.tail,log.p) = p
# however, there may be no such q, since we are limited to
# the range [lp,up) where
# lp = pco_sropt(lb,df1,df2,sropt,ope,lower.tail,log.p)
# up = pco_sropt(ub,df1,df2,sropt,ope,lower.tail,log.p)
# if p < lp we return lb;
# if q >= up, we return ub;
.qco_sropt <- function(p,df1,df2,z.s,ope,lower.tail=TRUE,log.p=FALSE,
												lb=0,ub=Inf) { # nocov start
	if ((lb > ub) || (is.infinite(lb)) || (min(lb,ub) < 0))
		stop("nonsensical lb and/or ub")

	eqv.p <- if (log.p) exp(p) else p

	if (eqv.p == 1)
		return(ifelse(lower.tail,Inf,0))
	if (eqv.p == 0)
		return(ifelse(lower.tail,0,Inf))
	if ((eqv.p < 0) || (eqv.p > 1))
		return (NaN)

	if (!missing(ope)) 
		z.s <- .deannualize(z.s,ope)

	# create a function increasing in its argument that
	# we wish to zero
	# do *not* pass on ope b/c this function is a tight loop
	if (lower.tail) {
		zerf <- function(q) {
			pco_sropt(q,df1=df1,df2=df2,z.s=z.s,lower.tail=lower.tail,log.p=log.p) - p
		}
	} else {
		zerf <- function(q) {
			p - pco_sropt(q,df1=df1,df2=df2,z.s=z.s,lower.tail=lower.tail,log.p=log.p)
		}
	}
	flb <- zerf(lb)
	if (flb > 0)
		return(lb)
	if (is.infinite(ub)) {
		ub <- 1 + lb
		fub <- zerf(ub)
		while ((fub < 0) && (!is.infinite(ub))) {
			ub <- 2 * ub
			fub <- zerf(ub)
		}
		if (is.infinite(ub) && (fub < 0)) 
			return(ub)
	} else {
		fub <- zerf(ub)
		if (fub < 0)
			return(ub)
	}

	ncp <- uniroot(zerf,interval=c(lb,ub),
								 f.lower=flb,f.upper=fub)
	retv <- ifelse(missing(ope),ncp$root,.annualize(ncp$root,ope))
	return(retv)
} # nocov end
#' @export 
qco_sropt <- Vectorize(.qco_sropt,
											vectorize.args = c("p","df1","df2","z.s"),
											SIMPLIFY = TRUE)
#UNFOLD

# co-F
# pco_f, qco_f#FOLDUP
#  ' @title The 'confidence distribution' for non-central F distribution.
#  '
#  ' @description 
#  '
#  ' Distribution function and quantile function for the 'confidence
#  ' distribution' of the non-central F distribution. This is just an inversion
#  ' to perform inference on \eqn{\lambda}{lambda} given observed statistic 
#  ' \eqn{F}.
#  '
#  ' @details
#  ' 
#  ' Suppose \eqn{x} follows a non-central F distribution
#  ' for known degrees of freedom, and 
#  ' unknown non-centrality parameter \eqn{\lambda}{lambda}. The 
#  ' 'confidence distribution' views \eqn{\lambda}{lambda} as a random
#  ' quantity once \eqn{x} is observed. As such, the CDF of
#  ' the confidence distribution is the same as that of the 
#  ' non-central F distribution (up to a flip of \code{lower.tail});
#  ' while the quantile function is used to compute confidence
#  ' intervals on \eqn{\lambda}{lambda} given \eqn{x}{x}.
#  '
#  ' @usage
#  '
#  ' pco_f(q,df1,df2,x,lower.tail=TRUE,log.p=FALSE) 
#  '
#  ' qco_f(p,df1,df2,x,lower.tail=TRUE,log.p=FALSE,lb=0,ub=Inf) 
#  '
#  ' @param q vector of quantiles.
#  ' @param p vector of probabilities.
#  ' @param x an observed statistic following a non-central F distribution.
#  ' @param log.p logical; if TRUE, probabilities p are given as \eqn{\mbox{log}(p)}{log(p)}.
#  ' @param lower.tail logical; if TRUE (default), probabilities are
#  '        \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}{P[X > x]}.
#  ' @param lb the lower bound for the output of \code{qco_f}.
#  ' @param ub the upper bound for the output of \code{qco_f}.
#  ' @inheritParams qf
#  ' @keywords distribution 
#  ' @return \code{pco_f} gives the distribution function, and
#  ' \code{qco_f} gives the quantile function.
#  '
#  ' Invalid arguments will result in return value \code{NaN} with a warning.
#  ' @aliases qco_f 
#  ' @seealso \code{\link{pf},\link{qf}}
#  ' @export 
#  ' @template etc
#  ' @note
#  ' When \code{lower.tail} is true, \code{pco_f} is monotonic increasing 
#  ' with respect to \code{q}, and decreasing in \code{x}; these are reversed
#  ' when \code{lower.tail} is false. Similarly, \code{qco_f} is increasing
#  ' in \code{sign(as.double(lower.tail) - 0.5) * p} and
#  ' \code{- sign(as.double(lower.tail) - 0.5) * x}.
#  '
#  ' @examples 
#  '
#  '
pco_f <- function(q,df1,df2,x,lower.tail=TRUE,log.p=FALSE) {
	# delegate
	retv <- pf(q=x,df1=df1,df2=df2,ncp=q,
						 lower.tail=!lower.tail,log.p=log.p)  # sic the tail reversal
	return(retv)
}
# create a scalar function that we later vectorize. 
# 
# this inverts pco_f; note that when lower.tail=TRUE,
# pco_f is increasing in q, but decreasing in x.
# pco_f only accepts non-negative q 
#
# here we try to find lb <= q < ub such that
# pco_f(q,df1,df2,x,ope,lower.tail,log.p) = p
# however, there may be no such q, since we are limited to
# the range [lp,up) where
# lp = pco_f(lb,df1,df2,x,ope,lower.tail,log.p)
# up = pco_f(ub,df1,df2,x,ope,lower.tail,log.p)
# if p < lp we return lb;
# if q >= up, we return ub;
.qco_f <- function(p,df1,df2,x,lower.tail=TRUE,log.p=FALSE,
												lb=0,ub=Inf) { # nocov start
	if ((lb > ub) || (is.infinite(lb)) || (min(lb,ub) < 0))
		stop("nonsensical lb and/or ub")

	eqv.p <- if (log.p) exp(p) else p

	if (eqv.p == 1)
		return(ifelse(lower.tail,Inf,0))
	if (eqv.p == 0)
		return(ifelse(lower.tail,0,Inf))
	if ((eqv.p < 0) || (eqv.p > 1))
		return (NaN)

	# create a function increasing in its argument that
	# we wish to zero
	# do *not* pass on ope b/c this function is a tight loop
	if (lower.tail) {
		zerf <- function(q) {
			pco_f(q,df1=df1,df2=df2,x=x,lower.tail=lower.tail,log.p=log.p) - p
		}
	} else {
		zerf <- function(q) {
			p - pco_f(q,df1=df1,df2=df2,x=x,lower.tail=lower.tail,log.p=log.p)
		}
	}
	flb <- zerf(lb)
	if (flb > 0)
		return(lb)
	if (is.infinite(ub)) {
		ub <- 1 + lb
		fub <- zerf(ub)
		while ((fub < 0) && (!is.infinite(ub))) {
			ub <- 2 * ub
			fub <- zerf(ub)
		}
		if (is.infinite(ub) && (fub < 0)) 
			return(ub)
	} else {
		fub <- zerf(ub)
		if (fub < 0)
			return(ub)
	}

	ncp <- uniroot(zerf,interval=c(lb,ub),
								 f.lower=flb,f.upper=fub)
	retv <- ncp$root
	return(retv)
} # nocov end
#  ' @export 
qco_f <- Vectorize(.qco_f,
									 vectorize.args = c("p","df1","df2","x"),
									 SIMPLIFY = TRUE)
#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
