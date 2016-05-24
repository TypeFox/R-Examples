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
#' @include distributions.r

# 2FIX:
# add britten-jones weights CI

########################################################################
# Tests 
########################################################################

# equality of SR#FOLDUP

# the guts of a SR equality test, hinging on different
# covariance functions.
.sr_eq_guts <- function(X,contrasts,vcov.func=vcov) {
	X <- na.omit(X)
	n <- dim(X)[1]
	p <- dim(X)[2]
	if (is.null(contrasts))
		contrasts <- .atoeplitz(c(1,-1,array(0,p-2)),c(1,array(0,p-2)))
	k <- dim(contrasts)[1]
	if (dim(contrasts)[2] != p)
		stop("size mismatch in 'X', 'contrasts'")

	retval <- sr_vcov(X,vcov.func=vcov.func)

	# store these:
	retval$n <- n
	retval$k <- k

	# the test statistic:
	retval$ESR <- contrasts %*% retval$SR
	retval$COC <- contrasts %*% retval$Ohat %*% t(contrasts)

	return(retval)
}
#' @title Paired test for equality of Sharpe ratio
#'
#' @description 
#'
#' Performs a hypothesis test of equality of Sharpe ratios of p assets
#' given paired observations.
#'
#' @details 
#'
#' Given \eqn{n} \emph{i.i.d.} observations of the excess returns of
#' \eqn{p} strategies, we test
#' \deqn{H_0: \frac{\mu_i}{\sigma_i} = \frac{\mu_j}{\sigma_j}, 1 \le i < j \le p}{H0: sr1 = sr2 = ...}
#' using the method of Wright, et. al. 
#' 
#' More generally, a matrix of constrasts, \eqn{E}{E} can be given, and we can
#' test
#' \deqn{H_0: E s = 0,}{H0: E s = 0,}
#' where \eqn{s}{s} is the vector of Sharpe ratios of the \eqn{p} strategies.
#' 
#' When \eqn{E}{E} consists of a single row (a single contrast), as is the
#' case when the default contrasts are used and only two strategies are
#' compared, then an approximate t-test can be performed against the
#' alternative hypothesis \eqn{H_a: E s > 0}{Ha: E s > 0}
#' 
#' Both chi-squared and F- approximations are supported; the former is
#' described by Wright. \emph{et. al.}, the latter by Leung and Wong.
#' 
#' @usage
#'
#' sr_equality_test(X,type=c("chisq","F","t"),
#'                  alternative=c("two.sided","less","greater"),
#'                  contrasts=NULL,
#'                  vcov.func=vcov) 
#'
#' @param X an \eqn{n \times p}{n x p} matrix of paired observations.
#' @param contrasts an \eqn{k \times p}{k x p} matrix of the contrasts
#         to test. This defaults to a matrix which tests sequential equality.
#' @param type which approximation to use. \code{"chisq"} is preferred when
#'        the returns are non-normal, but the approximation is asymptotic.
#'        the \code{"t"} test is only supported when \eqn{k = 1}{k = 1}.
#' @param alternative a character string specifying the alternative hypothesis,
#'        must be one of \code{"two.sided"} (default), \code{"greater"} or
#'        \code{"less"}. You can specify just the initial letter.
#'        This is only relevant for the \code{"t"} test.
#'        \code{"greater"} corresponds to \eqn{H_a: E s > 0}{Ha: E s > 0}.
#' @param vcov.func a function which takes a model of class lm (one of
#'        the form x ~ 1), and produces a variance-covariance matrix.
#'        The default is \code{\link{vcov}}, which produces a 'vanilla'
#'        estimate of covariance. Other sensible options are
#'        \code{vcovHAC} from the \code{sandwich} package.
#'
#' @keywords htest
#' @return Object of class \code{htest}, a list of the test statistic,
#' the size of \code{X}, and the \code{method} noted.
#' @seealso \code{\link{sr_test}}
#' @export 
#' @template etc
#' @template sr
#' @references 
#'
#' Wright, J. A., Yam, S. C. P., and Yung, S. P. "A note on the test for the
#' equality of multiple Sharpe ratios and its application on the evaluation
#' of iShares." J. Risk. to appear. 
#' \url{http://www.risk.net/journal-of-risk/technical-paper/2340067/a-test-for-the-equality-of-multiple-sharpe-ratios}
#'
#' Leung, P.-L., and Wong, W.-K. "On testing the equality of multiple Sharpe ratios, with 
#' application on the evaluation of iShares." J. Risk 10, no. 3 (2008): 15--30.
#' \url{http://papers.ssrn.com/sol3/papers.cfm?abstract_id=907270}
#'
#' Memmel, C. "Performance hypothesis testing with the Sharpe ratio." Finance
#' Letters 1 (2003): 21--23.
#'
#' @template ref-LW
#' @template ref-Lo
#'
#' @examples 
#' # under the null 
#' rv <- sr_equality_test(matrix(rnorm(500*5),ncol=5))
#' # under the alternative (but with identity covariance)
#' ope <- 253
#' nyr <- 10
#' nco <- 5
#' rets <- 0.01 * sapply(seq(0,1.7/sqrt(ope),length.out=nco),
#'   function(mu) { rnorm(ope*nyr,mean=mu,sd=1) })
#' rv <- sr_equality_test(rets)
#' 
#' \dontrun{
#' # using real data
#' if (require(quantmod)) {
#'   get.ret <- function(sym,...) {
#'     OHLCV <- getSymbols(sym,auto.assign=FALSE,...)
#'     lrets <- diff(log(OHLCV[,paste(c(sym,"Adjusted"),collapse=".",sep="")]))
#'     lrets[-1,]
#'   }
#'   get.rets <- function(syms,...) { some.rets <- do.call("cbind",lapply(syms,get.ret,...)) }
#'   some.rets <- get.rets(c("IBM","AAPL","NFLX","SPY"))
#'   pvs <- sr_equality_test(some.rets)
#' }
#' # test for uniformity
#' pvs <- replicate(1024,{ x <- sr_equality_test(matrix(rnorm(400*5),400,5),type="chisq")
#'                        x$p.value })
#' plot(ecdf(pvs))
#' abline(0,1,col='red') 
#' }
#' \dontrun{
#' if (require(sandwich)) {
#'   set.seed(as.integer(charToRaw("0b2fd4e9-3bdf-4e3e-9c75-25c6d18c331f")))
#'   n.manifest <- 10
#'   n.latent <- 4
#'   n.day <- 1024
#'   snr <- 0.95
#'   latent.rets <- matrix(rnorm(n.day*n.latent),ncol=n.latent) %*%
#'		matrix(runif(n.latent*n.manifest),ncol=n.manifest)
#'   noise.rets <- matrix(rnorm(n.day*n.manifest),ncol=n.manifest)
#'   some.rets <- snr * latent.rets + sqrt(1-snr^2) * noise.rets
#'   # naive vcov
#'   pvs0 <- sr_equality_test(some.rets)
#'   # HAC vcov
#'   pvs1 <- sr_equality_test(some.rets,vcov.func=vcovHAC)
#'   # more elaborately:
#'   pvs <- sr_equality_test(some.rets,vcov.func=function(amod) {
#'		vcovHAC(amod,prewhite=TRUE) })
#' }
#' }
#'
#'@export
sr_equality_test <- function(X,type=c("chisq","F","t"),
														 alternative=c("two.sided","less","greater"),
														 contrasts=NULL,
														 vcov.func=vcov) {
	# all this stolen from t.test.default:
	alternative <- match.arg(alternative)
	dname <- deparse(substitute(X))

	# delegate
	srmom <- .sr_eq_guts(X,contrasts,vcov.func=vcov.func)

	type <- match.arg(type)
	if ((type == "t") && (srmom$k != 1))
		stop("can only perform t-test on single contrast");

	if (type == "t") {
		#ts <- srmom$ESR * sqrt(srmom$n / srmom$COC)
		ts <- srmom$ESR * sqrt(1 / srmom$COC)
		names(ts) <- "t"
		pval <- switch(alternative,
									 two.sided = .oneside2two(pt(ts,df=srmom$n-1)),
									 less = pt(ts,df=srmom$n-1,lower.tail=TRUE),
									 greater = pt(ts,df=srmom$n-1,lower.tail=FALSE))
		statistic <- ts
	} else {
		#T2 <- srmom$n * t(srmom$ESR) %*% solve(srmom$COC,srmom$ESR)
		T2 <- t(srmom$ESR) %*% solve(srmom$COC,srmom$ESR)
		names(T2) <- "T2"
		pval <- switch(type,
									 chisq = pchisq(T2,df=srmom$k,ncp=0,lower.tail=FALSE),
									 F = pf((srmom$n-srmom$k) * T2/((srmom$n-1) * srmom$k),
													df1=srmom$k,df2=srmom$n-srmom$k,lower.tail=FALSE))
		statistic <- T2
		if (alternative != "two.sided") {
			warning("cannot perform directional tests on T^2")
			alternative <- "two.sided"
		}
	}

	# attach names
	names(srmom$k) <- "contrasts"
	method <- paste(c("test for equality of Sharpe ratio, via",type,"test"),collapse=" ")
	names(srmom$SR) <- sapply(1:srmom$p,function(x) { paste(c("strat",x),collapse="_") })

	cozeta <- 0
	names(cozeta) <- "sum squared contrasts of SNR"

	retval <- list(statistic = statistic, parameter = srmom$k,
							 df1 = srmom$p, df2 = srmom$n, p.value = pval, 
							 SR = srmom$SR, null.value = cozeta,
							 alternative = alternative,
							 method = method, data.name = dname)
	class(retval) <- "htest"
	return(retval)
}
#UNFOLD

# SR test#FOLDUP
#getAnywhere("t.test.default")
#getAnywhere("print.htest")
#' @title test for Sharpe ratio
#'
#' @description 
#'
#' Performs one and two sample tests of Sharpe ratio on vectors of data.
#'
#' @details 
#'
#' Given \eqn{n}{n} observations \eqn{x_i}{xi} from a normal random variable,
#' with mean \eqn{\mu}{mu} and standard deviation \eqn{\sigma}{sigma}, tests
#' \deqn{H_0: \frac{\mu}{\sigma} = S}{H0: mu/sigma = S}
#' against two or one sided alternatives.
#' 
#' Can also perform two sample tests of Sharpe ratio. For paired observations
#' \eqn{x_i}{xi} and \eqn{y_i}{yi}, tests
#' \deqn{H_0: \frac{\mu_x}{\sigma_x} = \frac{\mu_u}{\sigma_y}}{H0: mu_x sigma_y = mu_y sigma_x}
#' against two or one sided alternative, via 
#' \code{\link{sr_equality_test}}.
#'
#' For unpaired (and independent) observations, tests
#' \deqn{H_0: \frac{\mu_x}{\sigma_x} - \frac{\mu_u}{\sigma_y} = S}{H0: mu_x / sigma_x - mu_y / sigma_y = S}
#' against two or one-sided alternatives via the upsilon distribution.
#' 
#' @usage
#'
#' sr_test(x,y=NULL,alternative=c("two.sided","less","greater"),
#'         zeta=0,ope=1,paired=FALSE,conf.level=0.95)
#'
#' @param x a (non-empty) numeric vector of data values, or an
#'    object of class \code{sr}, containing a scalar sample Sharpe estimate.
#' @param y an optional (non-empty) numeric vector of data values, or
#'    an object of class \code{sr}, containing a scalar sample Sharpe estimate.
#'    Only an unpaired test can be performed when at least one of \code{x} and \code{y} are
#'    of class \code{sr}
#' @param alternative a character string specifying the alternative hypothesis,
#'       must be one of \code{"two.sided"} (default), \code{"greater"} or
#'       \code{"less"}.  You can specify just the initial letter.
#' @param zeta a number indicating the null hypothesis offset value, the
#' \eqn{S} value.
#' @template param-ope
#' @param paired a logical indicating whether you want a paired test.
#' @param conf.level confidence level of the interval. 
#' @template param-ellipsis
#' @keywords htest
#' @return A list with class \code{"htest"} containing the following components:
#' \item{statistic}{the value of the t- or Z-statistic.}
#' \item{parameter}{the degrees of freedom for the statistic.}
#' \item{p.value}{the p-value for the test.}
#' \item{conf.int}{a confidence interval appropriate to the specified alternative hypothesis. NYI for some cases.}
#' \item{estimate}{the estimated Sharpe or difference in Sharpes depending on whether it was a one-sample test or a two-sample test. Annualized}
#' \item{null.value}{the specified hypothesized value of the Sharpe or difference of Sharpes depending on whether it was a one-sample test or a two-sample test.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string indicating what type of test was performed.}
#' \item{data.name}{a character string giving the name(s) of the data.}
#' @seealso \code{\link{sr_equality_test}}, \code{\link{sr_unpaired_test}}, \code{\link{t.test}}.
#' @export 
#' @template etc
#' @template sr
#' @template ref-upsilon
#' @rdname sr_test
#' @examples 
#' # should reject null
#' x <- sr_test(rnorm(1000,mean=0.5,sd=0.1),zeta=2,ope=1,alternative="greater")
#' x <- sr_test(rnorm(1000,mean=0.5,sd=0.1),zeta=2,ope=1,alternative="two.sided")
#' # should not reject null
#' x <- sr_test(rnorm(1000,mean=0.5,sd=0.1),zeta=2,ope=1,alternative="less")
#'
#' # test for uniformity
#' pvs <- replicate(128,{ x <- sr_test(rnorm(1000),ope=253,alternative="two.sided")
#'                         x$p.value })
#' plot(ecdf(pvs))
#' abline(0,1,col='red') 
#' # testing an object of class sr
#' asr <- as.sr(rnorm(1000,1 / sqrt(253)),ope=253)
#' checkit <- sr_test(asr,zeta=0)
#'
#' @export
sr_test <- function(x,y=NULL,alternative=c("two.sided","less","greater"),
										zeta=0,ope=1,paired=FALSE,conf.level=0.95) {
	# much of this stolen from t.test.default:
	alternative <- match.arg(alternative)
	#if (is.sr(x) && is.null(y)) {
		#retv <- .sr_test_on_sr(z=x,alternative=alternative,
													 #zeta=zeta,conf.level=conf.level)
		#return(retv)
	#}
	# much of this stolen from t.test.default:
	if (!missing(zeta) && (length(zeta) != 1 || is.na(zeta))) 
		stop("'zeta' must be a single number")
	if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || 
		conf.level < 0 || conf.level > 1)) 
		stop("'conf.level' must be a single number between 0 and 1")

	if (!is.null(y)) {
		dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
		if (paired) {
			xok <- yok <- complete.cases(x, y)
		} else {
			yok <- !is.na(y)
			xok <- !is.na(x)
		}
		y <- y[yok]
		x <- x[xok]
	} else {
		dname <- deparse(substitute(x))
		if (paired) 
			warning("'y' is missing for paired test; assuming unpaired")
	}

	if (is.null(y)) {#FOLDUP
		# delegate
		if (is.sr(x)) { 
			srx <- x 
		} else {
			srx <- as.sr(x,c0=0,ope=ope,na.rm=TRUE)
		} 
		retv <- .sr_test_on_sr(srx,alternative=alternative,zeta=zeta,conf.level=conf.level)
		retv$data.name <- dname
		names(retv$estimate) <- paste(c("Sharpe ratio of ",dname),sep=" ",collapse="")
		return(retv)
	} #UNFOLD
	else {#FOLDUP
		ny <- length(y)
		if (paired) {#FOLDUP
			nx <- length(x)
			if (zeta != 0)
				stop("cannot test 'zeta' != 0 for paired test")
			if (nx != ny)
				stop("'x','y' must be same length")
			df <- nx - 1

			subtest <- sr_equality_test(cbind(x,y),type="t",alternative=alternative)
			# x minus y
			estimate <- - diff(as.vector(subtest$SR))
			estimate <- .annualize(estimate,ope)
			method <- "Paired sr-test"

			statistic <- subtest$statistic
			names(statistic) <- "t"
			# 2FIX: take df from subtest?
			pval <- subtest$p.value
		} #UNFOLD
		else {#FOLDUP
			srx <- as.sr(x,c0=0,ope=1,na.rm=TRUE)
			sry <- as.sr(y,c0=0,ope=1,na.rm=TRUE)
			retval <- sr_unpaired_test(list(srx,sry),c(1,-1),0,
																 alternative=alternative,
																 ope=ope,conf.level=conf.level)
			return(retval)
		}#UNFOLD
		names(estimate) <- "difference in Sharpe ratios"
	}#UNFOLD

	names(df) <- "df"
	names(zeta) <- "difference in signal-noise ratios"
	#attr(cint, "conf.level") <- conf.level
	retval <- list(statistic = statistic, parameter = df,
								 estimate = estimate, p.value = pval, 
								 alternative = alternative, null.value = zeta,
								 method = method, data.name = dname)
	class(retval) <- "htest"
	return(retval)
}
# screw it, this used to be sr_test.sr, but tired of fighting with roxygen
# and R CMD check --as-cran warnings.
.sr_test_on_sr <- function(z,alternative=c("two.sided","less","greater"),
											 zeta=0,conf.level=0.95) {
	# all this stolen from t.test.default:
	alternative <- match.arg(alternative)
	if (!missing(zeta) && (length(zeta) != 1 || is.na(zeta))) 
		stop("'zeta' must be a single number")
	if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || 
		conf.level < 0 || conf.level > 1)) 
		stop("'conf.level' must be a single number between 0 and 1")
	dname <- deparse(substitute(z))

	df <- z$df
	if (any(df < 1)) stop("not enough 'x' observations")
	statistic <- .sr2t(z)
	names(statistic) <- "t"
	estimate <- z$sr
	names(estimate) <- paste(c("Sharpe ratio of",dname),sep=" ",collapse="")

	method <- "One Sample sr test"

	if (alternative == "less") {
		pval <- .psr(z, zeta=zeta, lower.tail=TRUE)
		cint <- confint(z,type="exact",level.lo=0,level.hi=conf.level)
	}
	else if (alternative == "greater") {
		pval <- .psr(z, zeta=zeta, lower.tail=FALSE)
		cint <- confint(z,type="exact",level.lo=1-conf.level,level.hi=1)
	}
	else {
		pval <- .oneside2two(.psr(z, zeta=zeta, lower.tail=TRUE))
		cint <- confint(z,type="exact",level=conf.level)
	}

	names(df) <- "df"
	names(zeta) <- "signal-noise ratio"
	attr(cint, "conf.level") <- conf.level
	retval <- list(statistic = statistic, parameter = df,
								 estimate = estimate, p.value = pval, 
								 alternative = alternative, null.value = zeta,
								 method = method, data.name = dname)
	class(retval) <- "htest"
	return(retval)
}
#' @title test for equation on unpaired Sharpe ratios
#'
#' @description 
#'
#' Performs hypothesis tests on a single equation on k independent samples of Sharpe ratio.
#'
#' @details 
#'
#' For \eqn{1 \le j \le k}{1 <= j <= k}, suppose you have \eqn{n_j}{n_j}
#' observations of a normal random variable with mean \eqn{\mu_j}{mu_j} and
#' standard deviation \eqn{\sigma_j}{sigma_j}, with all observations
#' independent. Given constants \eqn{a_j}{a_j} and value \eqn{b}{b}, this
#' code tests the null hypothesis
#' \deqn{H_0: \sum_j a_j \frac{\mu_j}{\sigma_j} = b}{H0: sum_j a_j mu_j/sigma_j = b}
#' against two or one sided alternatives.
#' 
#' @usage
#'
#' sr_unpaired_test(srs,contrasts=NULL,null.value=0,
#'   alternative=c("two.sided","less","greater"),
#'   ope=NULL,conf.level=0.95)
#'
#' @param srs a (non-empty) list of objects of class \code{sr}, each containing
#'  a scalar sample Sharpe estimate. Or a single object of class \code{sr} with
#'  multiple Sharpe estimates. If the \code{sr} objects have different
#'  annualizations (\code{ope} parameters), a warning is thrown, since
#'  it is presumed that the contrasts all have the same units, but the
#'  test proceeds.
#' @param contrasts an array of the constrasts, the \eqn{a_j}{a_j} values.
#'  Defaults to \code{c(1,-1,1,...)}.
#' @param null.value the constant null value, the \eqn{b}{b}.
#'  Defaults to 0.
#' @param ope the number of observations per 'epoch'. For convenience of
#'   interpretation, The Sharpe ratio is typically quoted in 'annualized' 
#'   units for some epoch, that is, 'per square root epoch', though returns 
#'   are observed at a frequency of \code{ope} per epoch. 
#'   The default value is to take the same \code{ope} from the input \code{srs}
#'   object, if it is unambiguous. Otherwise, it defaults to 1, with a warning
#'   thrown.
#' @inheritParams sr_test
#' @template param-ellipsis
#' @keywords htest
#' @return A list with class \code{"htest"} containing the following components:
#' \item{statistic}{\code{NULL} here.}
#' \item{parameter}{a list of upsilon parameters.}
#' \item{p.value}{the p-value for the test.}
#' \item{conf.int}{a confidence interval appropriate to the specified alternative hypothesis.}
#' \item{estimate}{the estimated equation value, just the weighted sum of the sample Sharpe ratios. Annualized}
#' \item{null.value}{the specified hypothesized value of the sum of Sharpes.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string indicating what type of test was performed.}
#' \item{data.name}{a character string giving the name(s) of the data.}
#' @seealso \code{\link{sr_equality_test}}, \code{\link{sr_test}}, \code{\link{t.test}}.
#' @export 
#' @template etc
#' @template sr
#' @template ref-upsilon
#' @rdname sr_unpaired_test
#' @examples 
#' # basic usage
#' set.seed(as.integer(charToRaw("set the seed")))
#' # default contrast is 1,-1,1,-1,1,-1
#' etc <- sr_unpaired_test(as.sr(matrix(rnorm(1000*6,mean=0.02,sd=0.1),ncol=6)))
#' print(etc)
#'
#' etc <- sr_unpaired_test(as.sr(matrix(rnorm(1000*4,mean=0.0005,sd=0.01),ncol=4)),
#'   alternative='greater')
#' print(etc)
#'
#' etc <- sr_unpaired_test(as.sr(matrix(rnorm(1000*4,mean=0.0005,sd=0.01),ncol=4)),
#'   contrasts=c(1,1,1,1),null.value=-0.1,alternative='greater')
#' print(etc)
#'
#' inp <- list(as.sr(rnorm(500)),as.sr(runif(200)-0.5),
#'             as.sr(rnorm(30)),as.sr(rnorm(100)))
#' etc <- sr_unpaired_test(inp)
#'
#' inp <- list(as.sr(rnorm(500)),as.sr(rnorm(100,mean=0.2,sd=1)))
#' etc <- sr_unpaired_test(inp,contrasts=c(1,1),null.value=0.2)
#' etc$conf.int
#'
#' @export
sr_unpaired_test <- function(srs,contrasts=NULL,null.value=0,alternative=c("two.sided","less","greater"),
														 ope=NULL,conf.level=0.95) {
	# much of this stolen from t.test.default:
	alternative <- match.arg(alternative)
	if (!missing(null.value) && (length(null.value) != 1 || is.na(null.value))) 
		stop("'null.value' must be a single number")
	if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || 
		conf.level < 0 || conf.level > 1)) 
		stop("'conf.level' must be a single number between 0 and 1")
	dname <- deparse(substitute(srs))

	stopifnot(is.sr(srs) || (is.list(srs) && all(unlist(lapply(srs,is.sr)))))
	listof <- is.sr(srs[[1]])
	if (listof) {
		nterm <- length(srs)
		#should be a better way to do this:
		vals.sr <- unlist(lapply(srs,function(x) { x$sr }))
		vals.rescal <- unlist(lapply(srs,function(x) { x$rescal }))
		vals.df <- unlist(lapply(srs,function(x) { x$df }))
		vals.ope <- unlist(lapply(srs,function(x) { x$ope }))
	} else {
		nterm <- length(srs$sr)
		# 2FIX: possibly recycle these out to common length..
		vals.sr <- srs$sr
		vals.rescal <- srs$rescal
		vals.df <- srs$df 
		vals.ope <- srs$ope
	}
	if (is.null(contrasts)) { contrasts <- (-1)^(1 + seq_len(nterm)) }
	stopifnot(nterm == length(contrasts))

	# get rid of matrix, vector, etc:
	contrasts <- c(contrasts)
	vals.sr <- c(vals.sr)
	vals.ope <- c(vals.ope)

	# check ope
	uni.ope <- unique(vals.ope)
	if (length(uni.ope) != 1) {
		warning("ambiguous ope; check units of your contrasts!")
	}
	# inferr missing ope
	if (is.null(ope)) {
		if (length(uni.ope) == 1) {
			ope <- uni.ope
		} else {
			warning("ambiguous ope, guessing 1")
			ope <- 1.0
		}
	}

	# nb. the Sharpe is tstat * sqrt(ope) * rescal
	# deannualize
	vals.sr.de <- .deannualize(vals.sr,vals.ope)
	conval <- 1 / sqrt(sum((vals.rescal * contrasts) ^ 2))

	ups.t <- conval * (contrasts * vals.sr)
	ups.df <- c(vals.df)
	ups.val <- conval * null.value
	
	#2FIX: lower tail here?
	less.lt <- FALSE
	if (alternative == "less") {
		pval <- sadists::pupsilon(ups.val,df=ups.df,t=ups.t,lower.tail=less.lt)
		ciq <- c(0.0,conf.level)
	}
	else if (alternative == "greater") {
		pval <- sadists::pupsilon(ups.val,df=ups.df,t=ups.t,lower.tail=! less.lt)
		ciq <- c(1.0-conf.level,1.0)
	}
	else {
		pval <- .oneside2two(sadists::pupsilon(ups.val,df=ups.df,t=ups.t,lower.tail=less.lt))
		ciq <- 0.5 * (1 + conf.level * c(-1.0,1.0))
	}
	# figure out cis
	cint <- sadists::qupsilon(ciq,df=ups.df,t=ups.t)
	# convert back
	cint <- cint / conval
	cint <- .annualize(cint,ope)
	attr(cint, "conf.level") <- conf.level

	statistic <- NULL 
	#names(statistic) <- "null"
	names(null.value) <- "weighted sum of signal-noise ratios"

	estimate <- sum(contrasts * vals.sr)
	estimate <- .annualize(estimate,ope)
	method <- "unpaired k-sample sr-test"
	names(estimate) <- "equation on Sharpe ratios"

	# 2FIX:
	names(ups.df) <- "df"
	retval <- list(statistic = statistic, parameter = ups.df,
								 estimate = estimate, p.value = pval, 
								 alternative = alternative, null.value = null.value, conf.int = cint,
								 method = method, data.name = dname)
	class(retval) <- "htest"
	return(retval)
}
#' @title Power calculations for Sharpe ratio tests
#'
#' @description 
#'
#' Compute power of test, or determine parameters to obtain target power.
#'
#' @details 
#'
#' Suppose you perform a single-sample test for significance of the
#' Sharpe ratio based on the corresponding single-sample t-test. 
#' Given any three of: the effect size (the population SNR, \eqn{\zeta}{zeta}), 
#' the number of observations, and the type I and type II rates,
#' this function computes the fourth.
#'
#' This is a thin wrapper on \code{\link{power.t.test}}.
#'
#' Exactly one of the parameters \code{n}, \code{zeta}, \code{power}, and 
#' \code{sig.level} must be passed as NULL, and that parameter is determined 
#' from the others.  Notice that \code{sig.level} has non-NULL default, so NULL 
#' must be explicitly passed if you want to compute it.
#' 
#' @usage
#'
#' power.sr_test(n=NULL,zeta=NULL,sig.level=0.05,power=NULL,
#'               alternative=c("one.sided","two.sided"),ope=NULL) 
#'
#' @param n Number of observations
#' @param zeta the 'signal-to-noise' parameter, defined as the population
#'        mean divided by the population standard deviation, 'annualized'.
#' @param sig.level Significance level (Type I error probability).
#' @param power Power of test (1 minus Type II error probability).
#' @param alternative One- or two-sided test.
#' @template param-ope
#' @keywords htest
#' @return Object of class \code{power.htest}, a list of the arguments
#' (including the computed one) augmented with \code{method}, \code{note}
#' and \code{n.epoch} elements, the latter is the number of epochs 
#' under the given annualization (\code{ope}), \code{NA} if none given.
#' @seealso \code{\link{power.t.test}}, \code{\link{sr_test}}
#' @export 
#' @template etc
#' @template sr
#' @template ref-JW
#' @references 
#'
#' Lehr, R. "Sixteen S-squared over D-squared: A relation for crude 
#' sample size estimates." Statist. Med., 11, no 8 (1992): 1099--1102. 
#' doi: 10.1002/sim.4780110811
#'
#' @examples 
#' anex <- power.sr_test(253,1,0.05,NULL,ope=253) 
#' anex <- power.sr_test(n=253,zeta=NULL,sig.level=0.05,power=0.5,ope=253) 
#' anex <- power.sr_test(n=NULL,zeta=0.6,sig.level=0.05,power=0.5,ope=253) 
#' # Lehr's Rule 
#' zetas <- seq(0.1,2.5,length.out=51)
#' ssizes <- sapply(zetas,function(zed) { 
#'   x <- power.sr_test(n=NULL,zeta=zed,sig.level=0.05,power=0.8,
#'        alternative="two.sided",ope=253)
#'   x$n / 253})
#' # should be around 8.
#' print(summary(ssizes * zetas * zetas))
#' # e = n z^2 mnemonic approximate rule for 0.05 type I, 50% power
#' ssizes <- sapply(zetas,function(zed) { 
#'   x <- power.sr_test(n=NULL,zeta=zed,sig.level=0.05,power=0.5,ope=253)
#'   x$n / 253 })
#' print(summary(ssizes * zetas * zetas - exp(1)))
#' 
#'@export
power.sr_test <- function(n=NULL,zeta=NULL,sig.level=0.05,power=NULL,
													alternative=c("one.sided","two.sided"),
													ope=NULL) {
	# stolen from power.t.test
	if (sum(sapply(list(n, zeta, power, sig.level), is.null)) != 1) 
			stop("exactly one of 'n', 'zeta', 'power', and 'sig.level' must be NULL")
	if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > 
			sig.level | sig.level > 1)) 
			stop("'sig.level' must be numeric in [0, 1]")
	type <- "one.sample"
	alternative <- match.arg(alternative)
	if (!missing(ope) && !is.null(ope) && !is.null(zeta)) 
		zeta <- .deannualize(zeta,ope)
	# delegate
	# 2FIX: what happens if zeta = 0 ? bleah.
	subval <- power.t.test(n=n,delta=zeta,sd=1,sig.level=sig.level,
												 power=power,type=type,alternative=alternative,
												 strict=FALSE)
	# interpret
	subval$zeta <- subval$delta
	if (!missing(ope) && !is.null(ope)) {
		subval$zeta <- .annualize(subval$zeta,ope)
		subval$n.epoch <- subval$n / ope
	} else {
		subval$n.epoch <- NA
	}
	
	retval <- subval[c("n","n.epoch","zeta","sig.level","power","alternative","note","method")]
	retval <- structure(retval,class=class(subval))
	return(retval)
}
#UNFOLD

# SROPT test#FOLDUP
#' @title test for optimal Sharpe ratio
#'
#' @description 
#'
#' Performs one sample tests of Sharpe ratio of the Markowitz portfolio.
#'
#' @details 
#'
#' Suppose \eqn{x_i}{xi} are \eqn{n}{n} independent draws of a \eqn{q}{q}-variate
#' normal random variable with mean \eqn{\mu}{mu} and covariance matrix
#' \eqn{\Sigma}{Sigma}. This code tests the hypothesis
#' \deqn{H_0: \mu^{\top}\Sigma^{-1}\mu = \delta_0^2}{H0: mu' Sigma^-1 mu = delta_0^2}
#'
#' The default alternative hypothesis is the one-sided 
#' \deqn{H_1: \mu^{\top}\Sigma^{-1}\mu > \delta_0^2}{H1: mu' Sigma^-1 mu > delta_0^2}
#' but this can be set otherwise.
#' 
#' Note there is no 'drag' term here since this represents a linear offset of
#' the population parameter.
#' 
#' @usage
#'
#' sropt_test(X,alternative=c("greater","two.sided","less"),
#'             zeta.s=0,ope=1,conf.level=0.95)
#'
#' @param X a (non-empty) numeric matrix of data values, each row independent,
#'       each column representing an asset, or an object of 
#'       class \code{sropt}.
#' @param alternative a character string specifying the alternative hypothesis,
#'       must be one of \code{"two.sided"}, \code{"greater"} (default) or
#'       \code{"less"}.  You can specify just the initial letter.
#' @param zeta.s a number indicating the null hypothesis value.
#' @template param-ope
#' @param conf.level confidence level of the interval. (not used yet)
#' @keywords htest
#' @return A list with class \code{"htest"} containing the following components:
#' \item{statistic}{the value of the \eqn{T^2}-statistic.}
#' \item{parameter}{a list of the degrees of freedom for the statistic.}
#' \item{p.value}{the p-value for the test.}
#' \item{conf.int}{a confidence interval appropriate to the specified alternative hypothesis. NYI.}
#' \item{estimate}{the estimated optimal Sharpe, annualized}
#' \item{null.value}{the specified hypothesized value of the optimal Sharpe.}
#' \item{alternative}{a character string describing the alternative hypothesis.}
#' \item{method}{a character string indicating what type of test was performed.}
#' \item{data.name}{a character string giving the name(s) of the data.}
#' @seealso \code{\link{sr_test}}, \code{\link{t.test}}.
#' @export 
#' @template etc
#' @template sropt
#' @examples 
#'
#' # test for uniformity
#' pvs <- replicate(128,{ x <- sropt_test(matrix(rnorm(1000*4),ncol=4),alternative="two.sided")
#'                         x$p.value })
#' plot(ecdf(pvs))
#' abline(0,1,col='red') 
#' 
#' # input a sropt objects:
#' nfac <- 5
#' nyr <- 10
#' ope <- 253
#' # simulations with no covariance structure.
#' # under the null:
#' set.seed(as.integer(charToRaw("be determinstic")))
#' Returns <- matrix(rnorm(ope*nyr*nfac,mean=0,sd=0.0125),ncol=nfac)
#' asro <- as.sropt(Returns,drag=0,ope=ope)
#' stest <- sropt_test(asro,alternative="two.sided")
#'
#'@export
sropt_test <- function(X,alternative=c("greater","two.sided","less"),
										zeta.s=0,ope=1,conf.level=0.95) {
	# all this stolen from t.test.default:
	alternative <- match.arg(alternative)
	if (!missing(zeta.s) && (length(zeta.s) != 1 || is.na(zeta.s))) 
		stop("'zeta.s' must be a single number")
	if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || 
		conf.level < 0 || conf.level > 1)) 
		stop("'conf.level' must be a single number between 0 and 1")

	dname <- deparse(substitute(X))
	if (is.sropt(X)) {
		subtest <- X
		if (!is.null(ope))  
			subtest <- reannualize(subtest,new.ope=ope)
	} else {
		subtest <- as.sropt(X,ope=ope)
	}
	statistic <- subtest$T2
	names(statistic) <- "T2"
	estimate <- subtest$sropt
	names(estimate) <- "optimal Sharpe ratio of X"

	method <- "One Sample sropt test"

	df1 <- subtest$df1
	df2 <- subtest$df2

	# 2FIX: add CIs here.
	if (alternative == "less") {
		pval <- psropt(estimate, df1=df1, df2=df2, zeta.s=zeta.s, ope=ope)
	}
	else if (alternative == "greater") {
		pval <- psropt(estimate, df1=df1, df2=df2, zeta.s=zeta.s, ope=ope, lower.tail = FALSE)
	}
	else {
		pval <- .oneside2two(psropt(estimate, df1=df1, df2=df2, zeta.s=zeta.s, ope=ope))
	}

	names(df1) <- "df1"
	names(df2) <- "df2"
	df <- c(df1,df2)
	names(zeta.s) <- "optimal signal-noise ratio"
	#attr(cint, "conf.level") <- conf.level
	retval <- list(statistic = statistic, parameter = df,
								 estimate = estimate, p.value = pval, 
								 alternative = alternative, null.value = zeta.s,
								 method = method, data.name = dname)
	class(retval) <- "htest"
	return(retval)
}
#UNFOLD

# power of tests:#FOLDUP

## this could live on its own, but no longer belongs here.
#power.T2.test <- function(df1=NULL,df2=NULL,ncp=NULL,sig.level=0.05,power=NULL) {
	## stolen from power.anova.test
	#if (sum(sapply(list(df1, df2, ncp, power, sig.level), is.null)) != 1) 
		#stop("exactly one of 'df1', 'df2', 'ncp', 'power', and 'sig.level' must be NULL")
	#if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > 
		#sig.level | sig.level > 1)) 
		#stop("'sig.level' must be numeric in [0, 1]")
	#p.body <- quote({
		#delta2 <- df2 * ncp
		#pT2(qT2(sig.level, df1, df2, lower.tail = FALSE), df1, df2, delta2, lower.tail = FALSE)
	#})
	#if (is.null(power)) 
		#power <- eval(p.body)
	#else if (is.null(df1)) 
		#df1 <- uniroot(function(df1) eval(p.body) - power, c(1, 3e+03))$root
	#else if (is.null(df2)) 
		#df2 <- uniroot(function(df2) eval(p.body) - power, c(3, 1e+06))$root
	#else if (is.null(ncp))
		#ncp <- uniroot(function(ncp) eval(p.body) - power, c(0, 3e+01))$root
	#else if (is.null(sig.level)) 
		#sig.level <- uniroot(function(sig.level) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
	#else stop("internal error")
	#NOTE <- "one sided test"
	#METHOD <- "Hotelling test"
	#retval <- structure(list(df1 = df1, df2 = df2, ncp = ncp, delta2 = df2 * ncp,
								 #sig.level = sig.level, power = power, 
								 #note = NOTE, method = METHOD), class = "power.htest")
	#return(retval)
#}
#' @title Power calculations for optimal Sharpe ratio tests
#'
#' @description 
#'
#' Compute power of test, or determine parameters to obtain target power.
#'
#' @details 
#'
#' Suppose you perform a single-sample test for significance of the
#' optimal Sharpe ratio based on the corresponding single-sample T^2-test. 
#' Given any four of: the effect size (the population optimal SNR, 
#' \eqn{\zeta_*}{zeta*}), the number of assets, the number of observations, 
#' and the type I and type II rates, this function computes the fifth.
#'
#' Exactly one of the parameters \code{df1}, \code{df2}, 
#' \code{zeta.s}, \code{power}, and 
#' \code{sig.level} must be passed as NULL, and that parameter is determined 
#' from the others.  Notice that \code{sig.level} has non-NULL default, so NULL 
#' must be explicitly passed if you want to compute it.
#' 
#' @usage
#'
#' power.sropt_test(df1=NULL,df2=NULL,zeta.s=NULL,
#'                  sig.level=0.05,power=NULL,ope=1)
#'
#' @inheritParams dsropt
#' @param zeta.s the 'signal-to-noise' parameter, defined as ...
#' @param sig.level Significance level (Type I error probability).
#' @param power Power of test (1 minus Type II error probability).
#' @template param-ope
#' @keywords htest
#' @return Object of class \code{power.htest}, a list of the arguments
#' (including the computed one) augmented with \code{method}, \code{note}
#' and \code{n.epoch} elements, the latter is the number of epochs 
#' under the given annualization (\code{ope}), \code{NA} if none given.
#' @seealso \code{\link{power.t.test}}, \code{\link{sropt_test}}
#' @export 
#' @template etc
#' @template sropt
#'
#' @examples 
#' anex <- power.sropt_test(8,4*253,1,0.05,NULL,ope=253) 
#' 
#' @export
power.sropt_test <- function(df1=NULL,df2=NULL,zeta.s=NULL,
														 sig.level=0.05,power=NULL,ope=1) {
	# stolen from power.anova.test
	if (sum(sapply(list(df1, df2, zeta.s, power, sig.level), is.null)) != 1) 
		stop("exactly one of 'df1', 'df2', 'zeta.s', 'power', and 'sig.level' must be NULL")
	if (!is.null(sig.level) && !is.numeric(sig.level) || any(0 > 
		sig.level | sig.level > 1)) 
		stop("'sig.level' must be numeric in [0, 1]")
	p.body <- quote({
		psropt(qsropt(sig.level, df1, df2, ope, lower.tail = FALSE),
					 df1, df2, zeta.s, ope, lower.tail = FALSE)
	})
	if (is.null(power)) 
		power <- eval(p.body)
	else if (is.null(df1)) 
		df1 <- uniroot(function(df1) eval(p.body) - power, c(1, 3e+03))$root
	else if (is.null(df2)) 
		df2 <- uniroot(function(df2) eval(p.body) - power, c(3, 1e+06))$root
	else if (is.null(zeta.s))
		zeta.s <- uniroot(function(zeta.s) eval(p.body) - power, c(0, 3e+01))$root
	else if (is.null(sig.level)) 
		sig.level <- uniroot(function(sig.level) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
	else stop("internal error")
	NOTE <- "one sided test"
	METHOD <- "Hotelling test"
	retval <- structure(list(df1 = df1, df2 = df2, zeta.s = zeta.s, 
								 sig.level = sig.level, power = power, ope = ope, 
								 note = NOTE, method = METHOD), class = "power.htest")
	return(retval)
}

#UNFOLD

# spanning tests #FOLDUP
# ' @title spanning test for sub-portfolio.
# '
# ' @description 
# '
# ' Performs a test of portfolio spanning.
# '
# ' @details 
# '
# ' Suppose \eqn{x_i}{xi} are \eqn{n}{n} independent draws of a \eqn{q}{q}-variate
# ' normal random variable with mean \eqn{\mu}{mu} and covariance matrix
# ' \eqn{\Sigma}{Sigma}. Let \eqn{G} be some \eqn{n_g \times p}{n_g x p} matrix
# ' of rank \eqn{n_g}{n_g}. Let 
# ' \deqn{\Delta\zeta_{*} = \max_{w\,: G\Sigma w = 0}\frac{w^{\top}\mu}{\sqrt{w^{\top}\Sigma w}}}{Delta zeta* = max {(w'mu)/sqrt(w'Sigma w)|G Sigma w = 0}}
# '
# ' A spanning test tests the hypothesis
# ' \deqn{H_0: \Delta \zeta_{*} = 0}{H0: Delta zeta* = 0}
# ' against the alternative hypothesis
# ' \deqn{H_1: \Delta \zeta_{*} > 0}{H0: Delta zeta* > 0}
# '
# ' An alternative formulation, is as follows. Let
# ' \deqn{\zeta_{*,I} = \max_{w\,: w = I v}\frac{w^{\top}\mu}{\sqrt{w^{\top}\Sigma w}}}{zeta*,I = max {(w'mu)/sqrt(w'Sigma w)|w = I v}}
# ' and, similarly, let
# ' \deqn{\zeta_{*,G} = \max_{w\,: w = G v}\frac{w^{\top}\mu}{\sqrt{w^{\top}\Sigma w}}}{zeta*,G = max {(w'mu)/sqrt(w'Sigma w)|w = G v}}
# '
# ' Then
# ' \deqn{\left(\Delta\zeta_{*}\right)^2 = \zeta_{*,I}^2 - \zeta_{*,G}^2}{(Delta zeta*)^2 = zeta*,I^2 - zeta*,G^2}
# ' And so the spanning test tests 
# ' \deqn{H_0: \zeta_{*,I}^2 = \zeta_{*,G}^2}{H0: zeta*,I^2 = zeta*,G^2}
# ' against the alternative hypothesis
# ' \deqn{H_1: \zeta_{*,I}^2 > \zeta_{*,G}^2}{H0: zeta*,I^2 > zeta*,G^2}
# '
# ' The spanning test also provides estimates and confidence intervals on the
# ' difference 
# ' \eqn{\zeta_{*,I}^2 - \zeta_{*,G}^2}{zeta*,I^2 - zeta*,G^2}
# ' 
# ' @usage
# '
# ' # 2FIX: start here ... 
# '
# ' @template etc
# ' @template sropt
# ' @template ref-JW
# ' @export

#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
