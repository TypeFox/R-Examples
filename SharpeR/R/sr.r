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
# Created: 2013.04.16
# Copyright: Steven E. Pav, 2012-2013
# Author: Steven E. Pav
# Comments: Steven E. Pav

#' @include utils.r
#' @include distributions.r
#' @include estimation.r

########################################################################
# Sharpe Ratio#FOLDUP

# spawn a "SR" object.
#' @title Create an 'sr' object.
#'
#' @description 
#'
#' Spawns an object of class \code{sr}.
#'
#' @details
#'
#' The \code{sr} class contains information about a rescaled t-statistic.
#' The following are list attributes of the object:
#' \describe{
#' \item{sr}{The Sharpe ratio statistic.}
#' \item{df}{The d.f. of the equivalent t-statistic.}
#' \item{c0}{The drag 'risk free rate' used.}
#' \item{ope}{The 'observations per epoch'.}
#' \item{rescal}{The rescaling parameter.}
#' \item{epoch}{The string name of the 'epoch'.}
#' }
#'
#' The stored Sharpe statistic, \code{sr} is equal to the t-statistic 
#' times \eqn{rescal * sqrt{ope}}{rescal * sqrt(ope)}.
#'
#' For the most part, this constructor should \emph{not} be called directly,
#' rather \code{\link{as.sr}} should be called instead to compute the
#' Sharpe ratio.
#'
#' @usage
#'
#' sr(sr,df,c0=0,ope=1,rescal=sqrt(1/(df+1)),epoch="yr") 
#'
#' @param sr a Sharpe ratio statistic.
#' @param df the degrees of freedom of the equivalent t-statistic.
#' @param c0 the 'risk-free' or 'disastrous' rate of return. this is
#'        assumed to be given in the same units as x, \emph{not}
#'        in 'annualized' terms.
#' @template param-ope
#' @param rescal the rescaling parameter.
#' @param epoch the string representation of the 'epoch', defaulting
#'        to 'yr'.
#' @keywords univar 
#' @return a list cast to class \code{sr}.
#' @seealso \code{\link{as.sr}}
#' @rdname sr
#' @export sr
#' @template etc
#' @template sr
#'
#' @note
#' 2FIX: allow rownames? 
#'
#' @examples 
#' # roll your own.
#' ope <- 253
#' zeta <- 1.0
#' n <- 3 * ope
#' rvs <- rsr(1,n,zeta,ope=ope)
#' roll.own <- sr(sr=rvs,df=n-1,ope=ope,rescal=sqrt(1/n))
#' # put a bunch in. naming becomes a problem.
#' rvs <- rsr(5,n,zeta,ope=ope)
#' roll.own <- sr(sr=rvs,df=n-1,ope=ope,rescal=sqrt(1/n))
#'
sr <- function(sr,df,c0=0,ope=1,rescal=sqrt(1/(df+1)),epoch="yr") {
	retval <- list(sr = sr,df = df,c0 = c0,
								 ope = ope,rescal = rescal,epoch = epoch)
	class(retval) <- "sr"
	return(retval)
}
# compute SR in only one place. I hope.
.compute_sr <- function(mu,c0,sigma,ope) {
	sr <- (mu - c0) / sigma
	if (!missing(ope))
		sr <- sr * sqrt(ope)
	return(sr)
}
.get_strat_names <- function(xstat) {
	retval <- if (is.null(rownames(xstat))) rep(c("Sharpe"),length(xstat)) else rownames(xstat)
	return(retval)
}
#' @title Compute the Sharpe ratio.
#'
#' @description 
#'
#' Computes the Sharpe ratio of some observed returns.
#'
#' @details
#'
#' Suppose \eqn{x_i}{xi} are \eqn{n}{n} independent returns of some
#' asset.
#' Let \eqn{\bar{x}}{xbar} be the sample mean, and \eqn{s}{s} be
#' the sample standard deviation (using Bessel's correction). Let \eqn{c_0}{c0}
#' be the 'risk free rate'.  Then
#' \deqn{z = \frac{\bar{x} - c_0}{s}}{z = (xbar - c0)/s} 
#' is the (sample) Sharpe ratio.
#' 
#' The units of \eqn{z}{z} are \eqn{\mbox{time}^{-1/2}}{per root time}.
#' Typically the Sharpe ratio is \emph{annualized} by multiplying by
#' \eqn{\sqrt{\mbox{ope}}}{sqrt(ope)}, where \eqn{\mbox{ope}}{ope} 
#' is the number of observations
#' per year (or whatever the target annualization epoch.)
#'
#' Note that if \code{ope} is not given, the converter from \code{xts}
#' attempts to infer the observations per year, without regard to 
#' the name of the \code{epoch} given.
#'
#' @usage
#'
#' as.sr(x,c0=0,ope=1,na.rm=FALSE,epoch="yr")
#'
#' @param x vector of returns, or object of class \code{data.frame}, \code{xts},
#'        or \code{lm}.
#' @param c0 the 'risk-free' or 'disastrous' rate of return. this is
#'        assumed to be given in the same units as x, \emph{not}
#'        in 'annualized' terms.
#' @template param-ope
#' @param na.rm logical.  Should missing values be removed?
#' @param epoch the string representation of the 'epoch', defaulting
#'        to 'yr'.
#' @keywords univar 
#' @return a list containing the following components:
#' \describe{
#' \item{sr}{the annualized Sharpe ratio.}
#' \item{df}{the t-stat degrees of freedom.}
#' \item{c0}{the risk free term.}
#' \item{ope}{the annualization factor.}
#' \item{rescal}{the rescaling factor.}
#' \item{epoch}{the string epoch.}
#' }
#' cast to class \code{sr}.
#' @seealso sr-distribution functions, \code{\link{dsr}, \link{psr}, \link{qsr}, \link{rsr}}
#' @rdname as.sr
#' @export as.sr
#' @template etc
#' @template sr
#' @template ref-Lo
#'
#' @examples 
#' # Sharpe's 'model': just given a bunch of returns.
#' asr <- as.sr(rnorm(253*3),ope=253)
#' # or a matrix, with a name
#' my.returns <- matrix(rnorm(253*3),ncol=1)
#' colnames(my.returns) <- c("my strategy")
#' asr <- as.sr(my.returns)
#' # given an xts object:
#' \dontrun{
#' if (require(quantmod)) {
#'   IBM <- getSymbols('IBM',auto.assign=FALSE)
#'   lrets <- diff(log(IBM[,"IBM.Adjusted"]))
#'   asr <- as.sr(lrets,na.rm=TRUE)
#' }
#' }
#' # on a linear model, find the 'Sharpe' of the residual term
#' nfac <- 5
#' nyr <- 10
#' ope <- 253
#' set.seed(as.integer(charToRaw("determinstic")))
#' Factors <- matrix(rnorm(ope*nyr*nfac,mean=0,sd=0.0125),ncol=nfac)
#' Betas <- exp(0.1 * rnorm(dim(Factors)[2]))
#' Returns <- (Factors %*% Betas) + rnorm(dim(Factors)[1],mean=0.0005,sd=0.012)
#' APT_mod <- lm(Returns ~ Factors)
#' asr <- as.sr(APT_mod,ope=ope)
#' # try again, but make the Returns independent of the Factors.
#' Returns <- rnorm(dim(Factors)[1],mean=0.0005,sd=0.012)
#' APT_mod <- lm(Returns ~ Factors)
#' asr <- as.sr(APT_mod,ope=ope)
#'
#' # compute the Sharpe of a bunch of strategies:
#' my.returns <- matrix(rnorm(253*3*4),ncol=4)
#' asr <- as.sr(my.returns)  # without sensible colnames?
#' colnames(my.returns) <- c("strat a","strat b","strat c","strat d")
#' asr <- as.sr(my.returns)
#'   
as.sr <- function(x,c0=0,ope=1,na.rm=FALSE,epoch="yr") {
	UseMethod("as.sr", x)
}
.as.sr.unified <- function(x,mu,sigma,c0,ope,na.rm,epoch) {
	z <- .compute_sr(mu,c0,sigma,ope)
	dim(z) <- c(length(mu),1)

	# what bother; try to find the rownames
	if (length(dimnames(x)) == 1) 
		rnz <- unlist(dimnames(x))
	else 
		rnz <- unlist(dimnames(x)[2])
	if (is.null(rnz))
		rnz <- deparse(substitute(x))

	if (length(rnz) == dim(z)[1]) {
		rownames(z) <- rnz
	} else {
		rstub <- deparse(substitute(x))
		rownames(z) <- sapply(1:dim(z)[1],
													function(n) { sprintf('%s_%d',rstub,n) })
	}

	# make it a col vector if not otherwise defined.
	if (is.null(dim(x))) 
		dim(x) <- c(length(x),1)
	if (na.rm) 
		df <- apply(!is.na(x),2,sum)
	else
		df <- dim(x)[1]
	
	# n.b. sr$df is n-1
	retval <- sr(z,df=df-1,c0=c0,ope=ope,
							 rescal=1/sqrt(df),epoch=epoch)
	return(retval)
}
#' @rdname as.sr
#' @method as.sr default
#' @S3method as.sr default
as.sr.default <- function(x,c0=0,ope=1,na.rm=FALSE,epoch="yr") {
	mu <- mean(x,na.rm=na.rm)
	sigma <- sd(x,na.rm=na.rm)
	retval <- .as.sr.unified(x=x,mu=mu,sigma=sigma,c0=c0,ope=ope,
													 na.rm=na.rm,epoch=epoch)
	return(retval)
}
#' @rdname as.sr
#' @method as.sr matrix
#' @S3method as.sr matrix
as.sr.matrix <- function(x,c0=0,ope=1,na.rm=FALSE,epoch="yr") {
	mu <- apply(x,2,mean,na.rm=na.rm)
	sigma <- apply(x,2,sd,na.rm=na.rm)
	retval <- .as.sr.unified(x=x,mu=mu,sigma=sigma,c0=c0,ope=ope,
													 na.rm=na.rm,epoch=epoch)
	return(retval)
}
#' @rdname as.sr
#' @method as.sr data.frame
#' @S3method as.sr data.frame
as.sr.data.frame <- function(x,c0=0,ope=1,na.rm=FALSE,epoch="yr") {
	retval <- as.sr.matrix(x,c0=c0,ope=ope,na.rm=na.rm,epoch=epoch)
	return(retval)
}
#' @rdname as.sr
#' @method as.sr lm 
#' @S3method as.sr lm
as.sr.lm <- function(x,c0=0,ope=1,na.rm=FALSE,epoch="yr") {
	mu <- x$coefficients["(Intercept)"]
	sigma <- sqrt(deviance(x) / x$df.residual)
	z <- .compute_sr(mu,c0,sigma,ope)
	dim(z) <- c(1,1)
	rownames(z) <- deparse(substitute(x))
	XXinv <- vcov(x) / sigma^2
	rescal <- sqrt(XXinv["(Intercept)","(Intercept)"])
	retval <- sr(z,df=x$df.residual,c0=c0,ope=ope,
							 rescal=rescal,epoch=epoch)
	return(retval)
}
#' @rdname as.sr
#' @method as.sr xts 
#' @S3method as.sr xts
as.sr.xts <- function(x,c0=0,ope=1,na.rm=FALSE,epoch="yr") {
	if (missing(ope) && missing(epoch)) {
		ope <- .infer_ope_xts(x)
		epoch <- "yr"
	}
	retval <- as.sr.data.frame(as.data.frame(x),c0=c0,ope=ope,na.rm=na.rm,epoch=epoch)
	return(retval)
}
#' @rdname as.sr
#' @method as.sr timeSeries
#' @S3method as.sr timeSeries
as.sr.timeSeries <- function(x,c0=0,ope=1,na.rm=FALSE,epoch="yr") {
	# you want to do this, but requires xts package. oops.
	#retval <- as.sr.xts(as.xts(x),c0=c0,ope=ope,na.rm=na.rm,epoch=epoch)
	if (missing(ope) && missing(epoch)) {
		ope <- .infer_ope_xts(x)
		epoch <- "yr"
	}
	retval <- as.sr.data.frame(as.data.frame(x),c0=c0,ope=ope,na.rm=na.rm,epoch=epoch)
	return(retval)
}
#' @title Is this in the "sr" class?
#'
#' @description 
#'
#' Checks if an object is in the class \code{'sr'}
#'
#' @details
#'
#' To satisfy the minimum requirements of an S3 class.
#'
#' @usage
#'
#' is.sr(x)
#'
#' @param x an object of some kind.
#' @return a boolean.
#' @seealso sr
#' @template etc
#' @family sr
#' @export
#'
#' @examples 
#' rvs <- as.sr(rnorm(253*8),ope=253)
#' is.sr(rvs)
is.sr <- function(x) inherits(x,"sr")

# ' @S3method format sr
# ' @export
format.sr <- function(x,...) {
	# oh! ugly! ugly!
	retval <- capture.output(print.sr(x,...))
	return(retval)
}
#' @title Print values.
#'
#' @description 
#'
#' Displays an object, returning it \emph{invisibly}, 
#' (via \code{invisible(x)}.)
#'
#' @param x an object of class \code{sr} or \code{sropt}.
#' @template param-ellipsis
#'
#' @return the object, wrapped in \code{invisible}.
#' @rdname print
#' @method print sr
#' @S3method print sr
#' @export
#' @template etc
#' @template sr
#' @examples 
#' # compute a 'daily' Sharpe
#' mysr <- as.sr(rnorm(253*8),ope=1,epoch="day")
#' print(mysr)
#' # roll your own.
#' ope <- 253
#' zeta <- 1.0
#' n <- 6 * ope
#' rvs <- rsr(1,n,zeta,ope=ope)
#' roll.own <- sr(sr=rvs,df=n-1,ope=ope,rescal=sqrt(1/n))
#' print(roll.own)
#' # put a bunch in. naming becomes a problem.
#' rvs <- rsr(5,n,zeta,ope=ope)
#' roll.own <- sr(sr=rvs,df=n-1,ope=ope,rescal=sqrt(1/n))
#' print(roll.own)
#' # for sropt objects:
#' nfac <- 5
#' nyr <- 10
#' ope <- 253
#' # simulations with no covariance structure.
#' # under the null:
#' set.seed(as.integer(charToRaw("be determinstic")))
#' Returns <- matrix(rnorm(ope*nyr*nfac,mean=0,sd=0.0125),ncol=nfac)
#' asro <- as.sropt(Returns,drag=0,ope=ope)
#' print(asro)
print.sr <- function(x,...) {
	tval <- .sr2t(x)
	pval <- pt(tval,x$df,lower.tail=FALSE)
	serr <- se(x,type="t")
	coefs <- cbind(x$sr,serr,tval,pval)
	#colnames(coefs) <- c("stat","t.stat","p.value")
	colnames(coefs) <- c(paste(c("SR/sqrt(",x$epoch,")"),sep="",collapse=""),
											 "Std. Error","t value","Pr(>t)")
	rownames(coefs) <- .get_strat_names(x$sr)
	printCoefmat(coefs,P.values=TRUE,has.Pvalue=TRUE,
							 digits=max(2, getOption("digits") - 3),
							 cs.ind=c(1,2),tst.ind=c(3),dig.tst=2)
	# invisible return
	invisible(x)
}
# @hadley's suggested form
# print.sr <- function(x,...) cat(format(x,...), "\n")

# SR methods#FOLDUP
# 2FIX: make the x$rescal * sqrt(x$ope) occur in one place only ... 
# get the t-stat associated with an SR object.
.sr2t <- function(x) {
	tval <- x$sr / (x$rescal * sqrt(x$ope))
	return(tval)
}
# and the reverse
.t2sr <- function(x,tval) {
	srval <- tval * (x$rescal * sqrt(x$ope))
	return(srval)
}
.psr <- function(q,zeta,...) {
	retv <- prt(q$sr,df=q$df,K=(q$rescal * sqrt(q$ope)),rho=zeta,...)
	return(retv)
}
# no longer needed?
#.dsr <- function(q,zeta,...) {
	#retv <- drt(q$sr,df=q$df,K=(q$rescal * sqrt(q$ope)),rho=zeta,...)
	#return(retv)
#}

#' @title Change the annualization of a Sharpe ratio.
#'
#' @description 
#'
#' Changes the annualization factor of a Sharpe ratio statistic, or the rate at
#' which observations are made.
#'
#' @param object an object of class \code{sr} or \code{sropt}.
#' @param new.ope the new observations per epoch. If none given, it is
#' not updated.
#' @param new.epoch a string representation of the epoch. If none given, it is not
#' updated.
#' @return the input object with the annualization and/or epoch updated.
#' @seealso sr
#' @family sr
#' @seealso sropt
#' @family sropt
#' @template etc
#' @export
#' @rdname reannualize
#'
#' @examples 
#' # compute a 'daily' Sharpe
#' mysr <- as.sr(rnorm(253*8),ope=1,epoch="day")
#' # turn into annual 
#' mysr2 <- reannualize(mysr,new.ope=253,new.epoch="yr")
#'
#' # for sropt
#' ope <- 253
#' zeta.s <- 1.0  
#' df1 <- 10
#' df2 <- 6 * ope
#' rvs <- rsropt(1,df1,df2,zeta.s,ope,drag=0)
#' roll.own <- sropt(z.s=rvs,df1,df2,drag=0,ope=ope,epoch="yr")
#' # make 'monthly'
#' roll.monthly <- reannualize(roll.own,new.ope=21,new.epoch="mo.")
#' # make 'daily'
#' roll.daily <- reannualize(roll.own,new.ope=1,new.epoch="day")
reannualize <- function(object,new.ope=NULL,new.epoch=NULL) {
	UseMethod("reannualize", object)
}
#' @rdname reannualize
#' @method reannualize sr
#' @S3method reannualize sr
#' @export
reannualize.sr <- function(object,new.ope=NULL,new.epoch=NULL) {
	if (!is.sr(object)) stop("must give sr object")
	if (!missing(new.ope)) {
		object$sr <- object$sr * sqrt(new.ope / object$ope)
		object$ope <- new.ope
	}
	if (!missing(new.epoch)) object$epoch <- new.epoch
	return(object)
}
#' @rdname reannualize
#' @method reannualize sropt
#' @S3method reannualize sropt
#' @export
reannualize.sropt <- function(object,new.ope=NULL,new.epoch=NULL) {
	if (!is.sropt(object)) stop("must give sropt object")
	if (!missing(new.ope)) {
		object$sropt <- object$sropt * sqrt(new.ope / object$ope)
		object$ope <- new.ope
	}
	if (!missing(new.epoch)) object$epoch <- new.epoch
	return(object)
}
#UNFOLD
#UNFOLD

########################################################################
# Markowitz #FOLDUP

# markowitz, possibly on a subset, or basket, of the
# assets.
.sub.markowitz <- function(mu,Sigma,w=NULL,H=NULL) {
	if (!is.null(H)) {
		mu <- H %*% mu
		Sigma <- .qoform(H,Sigma)
	}
	if (is.null(w)) 
		w <- solve(Sigma,mu)
	zeta.sq <- t(mu) %*% w

	if (!is.null(H)) 
		w <- t(H) %*% mu
	retval <- list(w=w,zeta.sq=zeta.sq)
	return(retval)
}

markowitz <- function(mu,Sigma,df2,w=NULL,H=NULL) {
	# also computes Hotelling's T2
	retval <- .sub.markowitz(mu,Sigma,w,H)
	df1 <- ifelse(is.null(H),length(mu),dim(H)[1])
	retval <- c(retval,list(mu=mu,Sigma=Sigma,
						 df1=df1,df2=df2,
						 T2=df2*retval$zeta.sq))
	class(retval) <- c("markowitz")
	return(retval)
}

as.markowitz <- function(X,...) {
	UseMethod("as.markowitz", X)
}

# compute the (constrained) markowitz portfolio
as.markowitz.default <- function(X,mu=NULL,Sigma=NULL,...) {
	X <- na.omit(X)
	if (is.null(mu)) 
		mu <- colMeans(X)
	if (is.null(Sigma)) 
		Sigma <- cov(X)
	df2 <- dim(X)[1]
	retval <- markowitz(mu,Sigma,df2,...)
	return(retval)
}

#UNFOLD

########################################################################
# Optimal Sharpe ratio#FOLDUP

# spawn a "SROPT" object.
#' @title Create an 'sropt' object.
#'
#' @description 
#'
#' Spawns an object of class \code{sropt}.
#'
#' @details
#'
#' The \code{sropt} class contains information about a rescaled T^2-statistic.
#' The following are list attributes of the object:
#' \describe{
#' \item{sropt}{The (optimal) Sharpe ratio statistic.}
#' \item{df1}{The number of assets.}
#' \item{df2}{The number of observations.}
#' \item{drag}{The drag term, which is the 'risk free rate' divided by
#' the maximum risk.}
#' \item{ope}{The 'observations per epoch'.}
#' \item{epoch}{The string name of the 'epoch'.}
#' }
#'
#' For the most part, this constructor should \emph{not} be called directly,
#' rather \code{\link{as.sropt}} should be called instead to compute the
#' needed statistics.
#'
#' @usage
#'
#' sropt(z.s,df1,df2,drag=0,ope=1,epoch="yr",T2=NULL)
#'
#' @param z.s an optimum Sharpe ratio statistic.
#' @inheritParams dsropt
#' @param drag the 'drag' term, \eqn{c_0/R}{c0/R}. defaults to 0. It is assumed
#'        that \code{drag} has been annualized, \emph{i.e.} has been multiplied
#'        by \eqn{\sqrt{ope}}{sqrt(ope)}. This is in contrast to the \code{c0}
#'        term given to \code{\link{sr}}.
#' @param T2 the Hotelling \eqn{T^2}{T^2} statistic. If not given, it is
#'        computed from the given information.
#' @template param-ope
#' @param epoch the string representation of the 'epoch', defaulting
#'        to 'yr'.
#' @keywords univar 
#' @return a list cast to class \code{sropt}, with the following attributes:
#' \describe{
#' \item{sropt}{the optimal Sharpe statistic.}
#' \item{df1}{the number of assets.}
#' \item{df2}{the number of observed vectors.}
#' \item{drag}{the input \code{drag} term.}
#' \item{ope}{the input \code{ope} term.}
#' \item{epoch}{the input \code{epoch} term.}
#' \item{T2}{the Hotelling \eqn{T^2} statistic.}
#' }
#' @seealso \code{\link{as.sropt}}
#' @rdname sropt
#' @export 
#' @template etc
#' @template sropt
#'
#' @note
#' 2FIX: allow rownames?
#'
#' @examples 
#' # roll your own.
#' ope <- 253
#' zeta.s <- 1.0
#' df1 <- 10
#' df2 <- 6 * ope
#' set.seed(as.integer(charToRaw("fix seed")))
#' rvs <- rsropt(1,df1,df2,zeta.s,ope,drag=0)
#' roll.own <- sropt(z.s=rvs,df1,df2,drag=0,ope=ope)
#' print(roll.own)
#' # put a bunch in. naming becomes a problem.
#' rvs <- rsropt(5,df1,df2,zeta.s,ope,drag=0)
#' roll.own <- sropt(z.s=rvs,df1,df2,drag=0,ope=ope)
#' print(roll.own)
sropt <- function(z.s,df1,df2,drag=0,ope=1,epoch="yr",T2=NULL) {
	retval <- list(sropt = z.s,df1 = df1,df2 = df2,
								 drag = drag,ope = ope,epoch = epoch)
	if (is.null(T2)) 
		T2 <- .sropt2T(retval)
	retval$T2 <- T2
	# not sure I should do this ...
	retval$pval <- pT2(T2,retval$df1,retval$df2,lower.tail=FALSE)
	class(retval) <- "sropt"
	return(retval)
}
#' @title Compute the Sharpe ratio of the Markowitz portfolio.
#'
#' @description 
#'
#' Computes the Sharpe ratio of the Markowitz portfolio of some observed returns.
#'
#' @details
#' 
#' Suppose \eqn{x_i}{xi} are \eqn{n}{n} independent draws of a \eqn{q}{q}-variate
#' normal random variable with mean \eqn{\mu}{mu} and covariance matrix
#' \eqn{\Sigma}{Sigma}. Let \eqn{\bar{x}}{xbar} be the (vector) sample mean, and 
#' \eqn{S}{S} be the sample covariance matrix (using Bessel's correction). Let
#' \deqn{\zeta(w) = \frac{w^{\top}\bar{x} - c_0}{\sqrt{w^{\top}S w}}}{zeta(w) = (w'xbar - c0)/sqrt(w'Sw)}
#' be the (sample) Sharpe ratio of the portfolio \eqn{w}{w}, subject to 
#' risk free rate \eqn{c_0}{c0}.
#'
#' Let \eqn{w_*}{w*} be the solution to the portfolio optimization problem:
#' \deqn{\max_{w: 0 < w^{\top}S w \le R^2} \zeta(w),}{max {zeta(w) | 0 < w'Sw <= R^2},}
#' with maximum value \eqn{z_* = \zeta\left(w_*\right)}{z* = zeta(w*)}.
#' Then 
#' \deqn{w_* = R \frac{S^{-1}\bar{x}}{\sqrt{\bar{x}^{\top}S^{-1}\bar{x}}}}{%
#' w* = R S^-1 xbar / sqrt(xbar' S^-1 xbar)}
#' and
#' \deqn{z_* = \sqrt{\bar{x}^{\top} S^{-1} \bar{x}} - \frac{c_0}{R}}{%
#' z* = sqrt(xbar' S^-1 xbar) - c0/R}
#'
#' The units of \eqn{z_*}{z*} are \eqn{\mbox{time}^{-1/2}}{per root time}.
#' Typically the Sharpe ratio is \emph{annualized} by multiplying by
#' \eqn{\sqrt{\mbox{ope}}}{sqrt(ope)}, where \eqn{\mbox{ope}}{ope} 
#' is the number of observations
#' per year (or whatever the target annualization epoch.)
#'
#' Note that if \code{ope} and \code{epoch} are not given, the 
#' converter from \code{xts} attempts to infer the observations per epoch,
#' assuming yearly epoch.
#'
#' @usage
#'
#' as.sropt(X,drag=0,ope=1,epoch="yr")
#'
#' @param X matrix of returns, or \code{xts} object.
#' @inheritParams sropt 
#' @keywords univar 
#' @return An object of class \code{sropt}.
#' @seealso \code{\link{sropt}}, \code{\link{sr}}, sropt-distribution functions, 
#' \code{\link{dsropt}, \link{psropt}, \link{qsropt}, \link{rsropt}}
#' @rdname as.sropt
#' @export 
#' @template etc
#' @template sropt
#' @examples 
#' nfac <- 5
#' nyr <- 10
#' ope <- 253
#' # simulations with no covariance structure.
#' # under the null:
#' set.seed(as.integer(charToRaw("be determinstic")))
#' Returns <- matrix(rnorm(ope*nyr*nfac,mean=0,sd=0.0125),ncol=nfac)
#' asro <- as.sropt(Returns,drag=0,ope=ope)
#' # under the alternative:
#' Returns <- matrix(rnorm(ope*nyr*nfac,mean=0.0005,sd=0.0125),ncol=nfac)
#' asro <- as.sropt(Returns,drag=0,ope=ope)
#' # generating correlated multivariate normal data in a more sane way
#' if (require(MASS)) {
#'   nstok <- 10
#'   nfac <- 3
#'   nyr <- 10
#'   ope <- 253
#'   X.like <- 0.01 * matrix(rnorm(500*nfac),ncol=nfac) %*% 
#'     matrix(runif(nfac*nstok),ncol=nstok)
#'   Sigma <- cov(X.like) + diag(0.003,nstok)
#'   # under the null:
#'   Returns <- mvrnorm(ceiling(ope*nyr),mu=matrix(0,ncol=nstok),Sigma=Sigma)
#'   asro <- as.sropt(Returns,ope=ope)
#'   # under the alternative
#'   Returns <- mvrnorm(ceiling(ope*nyr),mu=matrix(0.001,ncol=nstok),Sigma=Sigma)
#'   asro <- as.sropt(Returns,ope=ope)
#' }
#' \dontrun{
#' # using real data.
#' if (require(quantmod)) {
#'   get.ret <- function(sym,...) {
#'     OHLCV <- getSymbols(sym,auto.assign=FALSE,...)
#'     lrets <- diff(log(OHLCV[,paste(c(sym,"Adjusted"),collapse=".",sep="")]))
#'     # chomp first NA!
#'     lrets[-1,]
#'   }
#'   get.rets <- function(syms,...) { some.rets <- do.call("cbind",lapply(syms,get.ret,...)) }
#'   some.rets <- get.rets(c("IBM","AAPL","A","C","SPY","XOM"))
#'   asro <- as.sropt(some.rets)
#' }
#' }
as.sropt <- function(X,drag=0,ope=1,epoch="yr") {
	UseMethod("as.sropt", X)
}
#' @rdname as.sropt
#' @method as.sropt default
#' @S3method as.sropt default
as.sropt.default <- function(X,drag=0,ope=1,epoch="yr") {
	# somehow call sropt!
	hotval <- as.markowitz(X)
	hotval <- c(hotval,list(ope=ope,drag=drag))
	# get the zeta.s
	z.s <- .T2sropt(hotval)
	dim(z.s) <- c(1,1)

	# this stinks
	retv <- sropt(z.s=z.s,df1=hotval$df1,df2=hotval$df2,
								drag=drag,ope=ope,epoch=epoch,T2=hotval$T2)

	# 2FIX: merge markowitz in?

	return(retv)
}
#' @rdname as.sropt
#' @method as.sropt xts
#' @S3method as.sropt xts
as.sropt.xts <- function(X,drag=0,ope=1,epoch="yr") {
	if (missing(ope)) {
		ope <- .infer_ope_xts(X)
	}
	retval <- as.sropt.default(X,drag=drag,ope=ope,epoch=epoch)
	return(retval)
}
#' @title Is this in the "sropt" class?
#'
#' @description 
#'
#' Checks if an object is in the class \code{'sropt'}
#'
#' @details
#'
#' To satisfy the minimum requirements of an S3 class.
#'
#' @usage
#'
#' is.sropt(x)
#'
#' @param x an object of some kind.
#' @return a boolean.
#' @seealso sropt
#' @template etc
#' @family sropt
#' @export
#'
#' @examples 
#' nfac <- 5
#' nyr <- 10
#' ope <- 253
#' # simulations with no covariance structure.
#' # under the null:
#' set.seed(as.integer(charToRaw("be determinstic")))
#' Returns <- matrix(rnorm(ope*nyr*nfac,mean=0,sd=0.0125),ncol=nfac)
#' asro <- as.sropt(Returns,drag=0,ope=ope)
#' is.sropt(asro)
is.sropt <- function(x) inherits(x,"sropt")

#' @rdname print
#' @method print sropt
#' @S3method print sropt
#' @export
print.sropt <- function(x,...) {
	Tval <- x$T2
	pval <- .sropt.pval(x)
	ci <- confint(x,level=0.95)
	coefs <- cbind(x$sropt,sric(x),ci,Tval,pval)
	colnames(coefs) <- c(paste(c("SR/sqrt(",x$epoch,")"),sep="",collapse=""),
											 paste(c("SRIC/sqrt(",x$epoch,")"),sep="",collapse=""),
											 colnames(ci)[1],colnames(ci)[2],
											 "T^2 value","Pr(>T^2)")
	rownames(coefs) <- .get_strat_names(x$sropt)
	printCoefmat(coefs,P.values=TRUE,has.Pvalue=TRUE,
							 digits=max(2, getOption("digits") - 3),
							 cs.ind=c(1,2,3,4),tst.ind=c(5),dig.tst=2)
}

# SROPT methods#FOLDUP
# get the T2-stat associated with an SROPT object.
.sropt2T <- function(x,Tval=x$T2) {
	if (is.null(Tval)) {
		tval <- .deannualize(x$sropt + x$drag,x$ope)
		Tval <- x$df2 * (tval^2)
	}
	return(Tval)
}
# and the reverse; with a sticky zero. ouch.
.T2sropt <- function(x,Tval=x$T2) {
	z.star <- sqrt(pmax(Tval,0) / x$df2)
	z.star <- .annualize(z.star,x$ope)
	z.star <- z.star - x$drag
	return(z.star)
}
.sropt.pval <- function(x,...) {
	Tval <- .sropt2T(x,...)
	retv <- x$pval
	if (is.null(retv)) 
		retv <- pT2(Tval,x$df1,x$df2,lower.tail=FALSE)
	return(retv)
}
#UNFOLD
#UNFOLD

########################################################################
# Delta Squared Optimal Sharpe ratio#FOLDUP

# spawn a "DEL_SROPT" object.
#' @title Create an 'del_sropt' object.
#'
#' @description 
#'
#' Spawns an object of class \code{del_sropt}.
#'
#' @details
#'
#' The \code{del_sropt} class contains information about the difference
#' between two rescaled T^2-statistics, useful for spanning
#' tests, and inference on hedged portfolios.
#' The following are list attributes of the object:
#' \describe{
#' \item{sropt}{The (optimal) Sharpe ratio statistic of
#' the 'full' set of assets.}
#' \item{sropt_sub}{The (optimal) Sharpe ratio statistic on
#' some subset, or linear subspace, of the assets.}
#' \item{df1}{The number of assets.}
#' \item{df2}{The number of observations.}
#' \item{df1.sub}{The number of degrees of freedom in the 
#' hedge constraint.}
#' \item{drag}{The drag term, which is the 'risk free rate' divided by
#' the maximum risk.}
#' \item{ope}{The 'observations per epoch'.}
#' \item{epoch}{The string name of the 'epoch'.}
#' }
#'
#' For the most part, this constructor should \emph{not} be called directly,
#' rather \code{\link{as.del_sropt}} should be called instead to compute the
#' needed statistics.
#'
#' @usage
#'
#' del_sropt(z.s,z.sub,df1,df2,df1.sub,drag=0,ope=1,epoch="yr")
#'
#' @param z.s an optimum Sharpe ratio statistic, on some set of assets.
#' @param z.sub an optimum Sharpe ratio statistic, on a linear subspace
#' of the assets.  If larger than \code{z.s} an error is thrown.
#' @param df1.sub the rank of the linear subspace of the hedge
#' constraint. 
#' by restricting attention to the subspace.
#' @inheritParams sropt
#' @template param-ope
#' @keywords univar 
#' @return a list cast to class \code{del_sropt}, with attributes
#' \describe{
#' \item{sropt}{the optimal Sharpe statistic.}
#' \item{sropt.sub}{the optimal Sharpe statistic on the subspace.}
#' \item{df1}{the number of assets.}
#' \item{df2}{the number of observed vectors.}
#' \item{df1.sub}{the input \code{df1.sub} term.}
#' \item{drag}{the input \code{drag} term.}
#' \item{ope}{the input \code{ope} term.}
#' \item{T2}{the Hotelling \eqn{T^2} statistic.}
#' \item{T2.sub}{the Hotelling \eqn{T^2} statistic on the subspace.}
#' }
#'
#' @seealso \code{\link{as.del_sropt}}
#' @rdname del_sropt
#' @export 
#' @template etc
#' @template del_sropt
#' @template warning
#'
#' @note
#' 2FIX: allow rownames?
#'
#' @examples 
#' # roll your own.
#' ope <- 253
#'
#' set.seed(as.integer(charToRaw("be determinstic")))
#' n.stock <- 10
#' X <- matrix(rnorm(1000*n.stock),nrow=1000)
#' Sigma <- cov(X)
#' mu <- colMeans(X)
#' w <- solve(Sigma,mu)
#' z <- t(mu) %*% w
#' n.sub <- 6
#' w.sub <- solve(Sigma[1:n.sub,1:n.sub],mu[1:n.sub])
#' z.sub <- t(mu[1:n.sub]) %*% w.sub
#' df1.sub <- n.stock - n.sub
#'
#' roll.own <- del_sropt(z.s=z,z.sub=z.sub,df1=10,df2=1000,
#'  df1.sub=df1.sub,ope=ope)
#' print(roll.own)
#'
del_sropt <- function(z.s,z.sub,df1,df2,df1.sub,drag=0,ope=1,epoch="yr") {
	retval <- list(sropt = z.s,sropt.sub = z.sub,
								 sropt.del = sqrt((z.s)^2 - (z.sub)^2),
								 df1 = df1,df2 = df2,df1.sub = df1.sub,
								 drag = drag,ope = ope,epoch = epoch)

	retval$T2.sub <- .sropt2T(list(sropt=z.sub,drag=drag,ope=ope,df2=df2))
	#retval$T2 <- .sropt2T(retval)
	retval$T2 <- .sropt2T(list(sropt=z.s,drag=drag,ope=ope,df2=df2))
	retval$T2.del <- retval$T2 - retval$T2.sub

	class(retval) <- "del_sropt"
	return(retval)
}
# the statistical guts.
.del_sropt.asF <- function(x) {
	# see Giri eqn (1.9) and section 3.
	# make Z ~ B(bdf1,bdf2) under the null
	nmin <- x$df2 - 1
	Z.num <- (1 + x$T2.sub/nmin)
	Z.denom <- (1 + x$T2/nmin)
	R1 <- 1 - (1 / Z.num)
	Z <- Z.num / Z.denom
	bdf1 <- (x$df2 - x$df1) / 2
	bdf2 <- (x$df1 - x$df1.sub) / 2
	# transform beta to F; 
	# define the F so that it is non-central under the alternative
	# (on the top)
	df1 <- 2*bdf2
	df2 <- 2*bdf1;
	Fval <- (df2 * (1-Z)) / (df1 * Z)
	pval <- pf(Fval,df1,df2,lower.tail=FALSE)
	# Giri's R1
	retval <- list(Fval=Fval,df1=df1,df2=df2,N=x$df2,pval=pval,R1=R1)
	return(retval)
}
#' @title Compute the Sharpe ratio of a hedged Markowitz portfolio.
#'
#' @description 
#'
#' Computes the Sharpe ratio of the hedged Markowitz portfolio of some observed returns.
#'
#' @details
#' 
#' Suppose \eqn{x_i}{xi} are \eqn{n}{n} independent draws of a \eqn{q}{q}-variate
#' normal random variable with mean \eqn{\mu}{mu} and covariance matrix
#' \eqn{\Sigma}{Sigma}. Let \eqn{G} be a \eqn{g \times q}{g x q} matrix
#' of rank \eqn{g}.
#' Let \eqn{\bar{x}}{xbar} be the (vector) sample mean, and 
#' \eqn{S}{S} be the sample covariance matrix (using Bessel's correction). 
#' Let
#' \deqn{\zeta(w) = \frac{w^{\top}\bar{x} - c_0}{\sqrt{w^{\top}S w}}}{zeta(w) = (w'xbar - c0)/sqrt(w'Sw)}
#' be the (sample) Sharpe ratio of the portfolio \eqn{w}{w}, subject to 
#' risk free rate \eqn{c_0}{c0}.
#'
#' Let \eqn{w_*}{w*} be the solution to the portfolio optimization 
#' problem:
#' \deqn{\max_{w: 0 < w^{\top}S w \le R^2,\,G S w = 0} \zeta(w),}{max {zeta(w) | 0 < w'Sw <= R^2, G S w = 0},}
#' with maximum value \eqn{z_* = \zeta\left(w_*\right)}{z* = zeta(w*)}.
#'
#' Note that if \code{ope} and \code{epoch} are not given, the 
#' converter from \code{xts} attempts to infer the observations per epoch,
#' assuming yearly epoch.
#'
#' @usage
#'
#' as.del_sropt(X,G,drag=0,ope=1,epoch="yr")
#'
#' @param X matrix of returns, or \code{xts} object.
#' @param G an \eqn{g \times q}{g x q} matrix of hedge constraints. A 
#' garden variety application would have \code{G} be one row of the
#' identity matrix, with a one in the column of the instrument to be
#' 'hedged out'.
#' @inheritParams sropt 
#' @keywords univar 
#' @return An object of class \code{del_sropt}.
#' @seealso \code{\link{del_sropt}}, \code{\link{sropt}}, 
#' \code{\link{sr}}
#' @rdname as.del_sropt
#' @export 
#' @template etc
#' @template del_sropt
#' @examples 
#' nfac <- 5
#' nyr <- 10
#' ope <- 253
#' # simulations with no covariance structure.
#' # under the null:
#' set.seed(as.integer(charToRaw("be determinstic")))
#' Returns <- matrix(rnorm(ope*nyr*nfac,mean=0,sd=0.0125),ncol=nfac)
#' # hedge out the first one:
#' G <- matrix(diag(nfac)[1,],nrow=1)
#' asro <- as.del_sropt(Returns,G,drag=0,ope=ope)
#' print(asro)
#' G <- diag(nfac)[c(1:3),]
#' asro <- as.del_sropt(Returns,G,drag=0,ope=ope)
#' # compare to sropt on the remaining assets
#' # they should be close, but not exact.
#' asro.alt <- as.sropt(Returns[,4:nfac],drag=0,ope=ope)
#' \dontrun{
#' # using real data.
#' if (require(quantmod)) {
#'   get.ret <- function(sym,...) {
#'     OHLCV <- getSymbols(sym,auto.assign=FALSE,...)
#'     lrets <- diff(log(OHLCV[,paste(c(sym,"Adjusted"),collapse=".",sep="")]))
#'     # chomp first NA!
#'     lrets[-1,]
#'   }
#'   get.rets <- function(syms,...) { some.rets <- do.call("cbind",lapply(syms,get.ret,...)) }
#'   some.rets <- get.rets(c("IBM","AAPL","A","C","SPY","XOM"))
#'   # hedge out SPY
#'   G <- diag(dim(some.rets)[2])[5,]
#'   asro <- as.del_sropt(some.rets,G)
#' }
#' }
as.del_sropt <- function(X,G,drag=0,ope=1,epoch="yr") {
	UseMethod("as.del_sropt",X)
}
#' @rdname as.del_sropt
#' @method as.del_sropt default
#' @S3method as.del_sropt default
as.del_sropt.default <- function(X,G,drag=0,ope=1,epoch="yr") {
	# somehow call sropt!
	hotval <- as.markowitz(X)
	hotval <- c(hotval,list(ope=ope,drag=drag))
	# get the zeta.s
	z.s <- .T2sropt(hotval)
	dim(z.s) <- c(1,1)
	if (!is.matrix(G)) G <- matrix(G,nrow=1)

	if (length(hotval$w) != dim(G)[2]) stop("wrong size hedge constraint?")

	hotval.sub <- markowitz(hotval$mu,hotval$Sigma,hotval$df2,H=G)
	hotval.sub <- c(hotval.sub,list(ope=ope,drag=drag))
	# get the zeta.sub
	z.sub <- .T2sropt(hotval.sub)
	dim(z.sub) <- c(1,1)

	# this stinks
	retv <- del_sropt(z.s=z.s,z.sub=z.sub,df1=hotval$df1,df2=hotval$df2,
										df1.sub=dim(G)[1],drag=drag,ope=ope,epoch=epoch)

	return(retv)
}
#' @rdname as.del_sropt
#' @method as.del_sropt xts
#' @S3method as.del_sropt xts
as.del_sropt.xts <- function(X,G,drag=0,ope=1,epoch="yr") {
	if (missing(ope)) {
		ope <- .infer_ope_xts(X)
	}
	retval <- as.del_sropt.default(X,G,drag=drag,ope=ope,epoch=epoch)
	return(retval)
}
#' @title Is this in the "del_sropt" class?
#'
#' @description 
#'
#' Checks if an object is in the class \code{'del_sropt'}
#'
#' @details
#'
#' To satisfy the minimum requirements of an S3 class.
#'
#' @usage
#'
#' is.del_sropt(x)
#'
#' @param x an object of some kind.
#' @return a boolean.
#' @seealso del_sropt
#' @template etc
#' @family del_sropt
#' @export
#'
#' @examples 
#' roll.own <- del_sropt(z.s=2,z.sub=1,df1=10,df2=1000,df1.sub=3,ope=1,epoch="yr")
#' is.sropt(roll.own)
is.del_sropt <- function(x) inherits(x,"del_sropt")

#' @rdname print
#' @method print del_sropt
#' @S3method print del_sropt
#' @export
print.del_sropt <- function(x,...) {
	Fandp <- .del_sropt.asF(x)
	ci <- confint(x,level=0.95)
	coefs <- cbind(x$sropt.del,ci,Fandp$Fval,Fandp$pval)
	colnames(coefs) <- c(paste(c("SR/sqrt(",x$epoch,")"),sep="",collapse=""),
											colnames(ci)[1],colnames(ci)[2],
											 "F value","Pr(>F)")
	rownames(coefs) <- .get_strat_names(x$sropt)
	printCoefmat(coefs,P.values=TRUE,has.Pvalue=TRUE,
							 digits=max(2, getOption("digits") - 3),
							 cs.ind=c(1),tst.ind=c(2),dig.tst=2)
}
#UNFOLD

## test them:
#X <- matrix(rnorm(1000*10,0.1),ncol=10)
#G <- diag(dim(X)[2])[1,]
#foo <- as.del_sropt(X,G)
#print(foo)

#X <- matrix(rnorm(1000*10,0.00001),ncol=10)
#foo <- as.del_sropt(X,G)
#print(foo)

# create a summary#FOLDUP
#' @title Summarize a Sharpe, or (delta) optimal Sharpe object.
#'
#' @description 
#'
#' Computes a \sQuote{summary} of an object, adding in some statistics.
#'
#' @details
#'
#' Enhances an object of class \code{sr}, \code{sropt} or \code{del_sropt} to also 
#' include t- or T-statistics, p-values, and so on.
#' 
#' @usage
#'
#' summary(obj)
#'
#' @param obj an object of class \code{sr}, \code{sropt} or \code{del_sropt}.
#' @return When an \code{sr} object is input, the object cast to class \code{summary.sr} with some
#' additional fields:
#' \describe{
#' \item{tval}{the equivalent t-statistic.}
#' \item{pval}{the p-value under the null.}
#' \item{serr}{the standard error of the Sharpe ratio.}
#' }
#' When an \code{sropt} object is input, the object cast to class \code{summary.sropt} with some
#' additional fields:
#' \describe{
#' \item{pval}{the p-value under the null.}
#' \item{SRIC}{the SRIC value, see \code{\link{sric}}.}
#' }
#'
#' @seealso \code{\link{print.sr}}.
#' @rdname summary
#' @export summary
#' @template etc
#' @template sr
#'
#' @examples 
#' # Sharpe's 'model': just given a bunch of returns.
#' set.seed(1234)
#' asr <- as.sr(rnorm(253*3),ope=253)
#' summary(asr)
summary <- function(obj) {
	UseMethod("summary", obj)
}
#' @rdname summary
#' @method summary sr
#' @S3method summary sr
summary.sr <- function(obj) {
	obj$tval <- .sr2t(obj)
	obj$pval <- pt(obj$tval,obj$df,lower.tail=FALSE)
	obj$serr <- se(obj,type="t")
	class(obj) <- "summary.sr"
	obj
}
#' @rdname summary
#' @method summary sropt
#' @S3method summary sropt
summary.sropt <- function(obj) {
	obj$pval <- .sropt.pval(obj)
	obj$SRIC <- sric(obj)
	#... anything else? KRS?
	class(obj) <- "summary.sropt"
	obj
}

#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
