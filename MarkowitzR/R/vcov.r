# Copyright 2014-2014 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav
#
# This file is part of MarkowitzR.
#
# MarkowitzR is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MarkowitzR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MarkowitzR.  If not, see <http://www.gnu.org/licenses/>.

# Created: 2014.01.31
# Copyright: Steven E. Pav, 2014
# Author: Steven E. Pav
# Comments: Steven E. Pav

#source("utils.r")

#' @title Compute variance covariance of 'Unified' Second Moment 
#'
#' @description
#'
#' Computes the variance covariance matrix of sample mean and second moment.
#'
#' @details
#'
#' Given \eqn{p}-vector \eqn{x}, the 'unified' sample is the 
#' \eqn{(p+1)(p+2)/2} vector of 1, \eqn{x}, and 
#' \eqn{\mbox{vech}(x x^{\top})}{vech(x x')} stacked on top
#' of each other.
#' Given \eqn{n} contemporaneous observations of \eqn{p}-vectors,
#' stacked as rows in the \eqn{n \times p}{n x p} matrix \eqn{X},
#' this function computes the mean and the variance-covariance
#' matrix of the 'unified' sample. 
#'
#' One may use the default method for computing covariance,
#' via the \code{\link{vcov}} function, or via a 'fancy' estimator,
#' like \code{sandwich:vcovHAC}, \code{sandwich:vcovHC}, \emph{etc.}
#'
#' @usage
#'
#' theta_vcov(X,vcov.func=vcov,fit.intercept=TRUE)
#'
#' @param X an \eqn{n \times p}{n x p} matrix of observed returns.
#' @param vcov.func a function which takes an object of class \code{lm},
#' and computes a variance-covariance matrix. If equal to the string
#' \code{"normal"}, we assume multivariate normal returns.
#' @param fit.intercept a boolean controlling whether we add a column
#' of ones to the data, or fit the raw uncentered second moment.
#' For now, must be true when assuming normal returns.
#' @keywords univar 
#'
#' @return a list containing the following components:
#' \item{mu}{a \eqn{q = (p+1)(p+2)/2} vector of 1, then the mean, 
#' then the vech'd second moment of the sample data.}
#' \item{Ohat}{the \eqn{q \times q}{q x q} estimated variance covariance 
#' matrix. When \code{fit.intercept} is true, the left column and top row
#' are all zeros.}
#' \item{n}{the number of rows in \code{X}.}
#' \item{pp}{the number of assets plus \code{as.numeric(fit.intercept)}.}
#' @seealso \code{\link{itheta_vcov}}.
#' @rdname theta_vcov
#' @export 
#' @template etc
#' @template ref-SEP13
#'
#' @note
#'
#' Replaces similar functionality from SharpeR package, but with 
#' modified API.
#'
#' @examples 
#' X <- matrix(rnorm(1000*3),ncol=3)
#' Sigmas <- theta_vcov(X)
#' Sigmas.n <- theta_vcov(X,vcov.func="normal")
#' Sigmas.n <- theta_vcov(X,fit.intercept=FALSE)
#'
#' # make it fat tailed:
#' X <- matrix(rt(1000*3,df=5),ncol=3)
#' Sigmas <- theta_vcov(X)
#' \dontrun{
#' if (require(sandwich)) {
#'  Sigmas <- theta_vcov(X,vcov.func=vcovHC)
#' }
#' }
#' # add some autocorrelation to X
#' Xf <- filter(X,c(0.2),"recursive")
#' colnames(Xf) <- colnames(X)
#' Sigmas <- theta_vcov(Xf)
#' \dontrun{
#' if (require(sandwich)) {
#'	Sigmas <- theta_vcov(Xf,vcov.func=vcovHAC)
#' }
#' }
#'
theta_vcov <- function(X,vcov.func=vcov,fit.intercept=TRUE) {
	set.coln(X)
	X <- na.omit(X)
	if (fit.intercept) 
		X <- cbind(1,X)
	p <- dim(X)[2]
	TF <- outer(1:p,1:p,FUN=">=")
	Y <- X[,col(TF)[TF]] * X[,row(TF)[TF]]
	Ynames <- paste(colnames(X)[col(TF)[TF]],"x",colnames(X)[row(TF)[TF]])
	colnames(Y) <- Ynames

	# model Y and estimate Sigma hat#FOLDUP
	if (is.character(vcov.func)) {
		n <- dim(X)[1]
		if (vcov.func != "normal") 
			stop("only understand function or 'normal'")
		if (!fit.intercept) 
			stop("unknown form for gaussian zero mean case?")

		mu <- colMeans(Y)
		rm(Y)
		Theta <- ivech(mu)

		iTheta <- solve(Theta)
		Dupp <- matrixcalc::duplication.matrix(p)
		Elim <- matrixcalc::elimination.matrix(p)
		# 2FIX there's a faster way to perform elimination, BTW!
		LXD <- Elim %*% (iTheta %x% iTheta) %*% Dupp
		DLXDU <- Dupp %*% LXD[,2:dim(LXD)[2]]
		
		inner <- (Theta %x% Theta) 
		H <- t(DLXDU) %*% inner %*% DLXDU
		Ohat <- (2 / n) * solve(H)
		# now gate with zeros:
		Ohat <- cbind(0,rbind(0,Ohat))
	} else {
		# take first moment of Y, then get standard error
		mod2 <- lm(Y ~ 1)
		rm(Y)
		mu <- mod2$coefficients
		Ohat <- vcov.func(mod2)
		rm(mod2)
	}#UNFOLD

	# make nice output
	dim(mu) <- c(length(mu),1)
	rownames(mu) <- Ynames
	dimnames(Ohat) <- list(Ynames,Ynames)
	retval <- list(mu=mu,Ohat=Ohat,n=dim(X)[1],pp=p)
	return(retval)
}

#' @title Compute variance covariance of Inverse 'Unified' Second Moment 
#'
#' @description
#'
#' Computes the variance covariance matrix of the inverse unified 
#' second moment matrix.
#'
#' @details
#'
#' Given \eqn{p}-vector \eqn{x} with mean \eqn{\mu}{mu} and
#' covariance, \eqn{\Sigma}{Sigma}, let \eqn{y} be \eqn{x}
#' with a one prepended. Then let 
#' \eqn{\Theta = E\left(y y^{\top}\right)}{Theta = E[yy']},
#' the uncentered second moment matrix. The inverse of
#' \eqn{\Theta}{Theta} contains the (negative) Markowitz portfolio 
#' and the precision matrix. 
#'
#' Given \eqn{n} contemporaneous observations of \eqn{p}-vectors,
#' stacked as rows in the \eqn{n \times p}{n x p} matrix \eqn{X},
#' this function estimates the mean and the asymptotic 
#' variance-covariance matrix of \eqn{\Theta^{-1}}{Theta^-1}.
#'
#' One may use the default method for computing covariance,
#' via the \code{\link{vcov}} function, or via a 'fancy' estimator,
#' like \code{sandwich:vcovHAC}, \code{sandwich:vcovHC}, \emph{etc.}
#'
#' @usage
#'
#' itheta_vcov(X,vcov.func=vcov,fit.intercept=TRUE)
#'
#'
#' @inheritParams theta_vcov
#' @return a list containing the following components:
#' \item{mu}{a \eqn{q = (p+1)(p+2)/2} vector of 1 + squared maximum
#' Sharpe, the negative Markowitz 
#' portfolio, then the vech'd precision matrix of the sample data}
#' \item{Ohat}{the \eqn{q \times q}{q x q} estimated variance 
#' covariance matrix.}
#' \item{n}{the number of rows in \code{X}.}
#' \item{pp}{the number of assets plus \code{as.numeric(fit.intercept)}.}
#' @keywords univar 
#' @seealso \code{\link{theta_vcov}}, \code{\link{mp_vcov}}.
#' @rdname itheta_vcov
#' @export 
#' @template etc
#' @template ref-SEP13
#'
#' @note
#'
#' By flipping the sign of \eqn{X}, the inverse of 
#' \eqn{\Theta}{Theta} contains the \emph{positive} Markowitz
#' portfolio and the precision matrix on \eqn{X}. Performing
#' this transform before passing the data to this function
#' should be considered idiomatic.
#'
#' @note
#'
#' A more general form of this function exists as \code{\link{mp_vcov}}.
#'
#' @note
#'
#' Replaces similar functionality from SharpeR package, but with 
#' modified API.
#'
#' @examples 
#' X <- matrix(rnorm(1000*3),ncol=3)
#' # putting in -X is idiomatic:
#' ism <- itheta_vcov(-X)
#' iSigmas.n <- itheta_vcov(-X,vcov.func="normal")
#' iSigmas.n <- itheta_vcov(-X,fit.intercept=FALSE)
#' # compute the marginal Wald test statistics:
#' qidx <- 2:ism$pp
#' wald.stats <- ism$mu[qidx] / sqrt(diag(ism$Ohat[qidx,qidx]))
#'
#' # make it fat tailed:
#' X <- matrix(rt(1000*3,df=5),ncol=3)
#' ism <- itheta_vcov(X)
#' qidx <- 2:ism$pp
#' wald.stats <- ism$mu[qidx] / sqrt(diag(ism$Ohat[qidx,qidx]))
#' \dontrun{
#' if (require(sandwich)) {
#'  ism <- itheta_vcov(X,vcov.func=vcovHC)
#'  qidx <- 2:ism$pp
#'  wald.stats <- ism$mu[qidx] / sqrt(diag(ism$Ohat[qidx,qidx]))
#' }
#' }
#' # add some autocorrelation to X
#' Xf <- filter(X,c(0.2),"recursive")
#' colnames(Xf) <- colnames(X)
#' ism <- itheta_vcov(Xf)
#' qidx <- 2:ism$pp
#' wald.stats <- ism$mu[qidx] / sqrt(diag(ism$Ohat[qidx,qidx]))
#' \dontrun{
#' if (require(sandwich)) {
#'	ism <- itheta_vcov(Xf,vcov.func=vcovHAC)
#'  qidx <- 2:ism$pp
#'  wald.stats <- ism$mu[qidx] / sqrt(diag(ism$Ohat[qidx,qidx]))
#' }
#' }
itheta_vcov <- function(X,vcov.func=vcov,fit.intercept=TRUE) {
	# delegate
	sm_est <- theta_vcov(X,vcov.func,fit.intercept)

	# interpret
	pp <- sm_est$pp
	Theta <- ivech(sm_est$mu)
	iTheta <- solve(Theta)
	iiTheta <- iTheta %x% iTheta  # kron 
	
	# elimination & duplication
	elim.idx <- lower.tri(Theta,diag=TRUE)
	dim(elim.idx) <- c(length(elim.idx),1)
	Dupp <- matrixcalc::duplication.matrix(pp)

	# delta method outer deriv
	H.mat <- iiTheta[elim.idx,] %*% Dupp
	Ohat <- (H.mat) %*% sm_est$Ohat %*% t(H.mat)

	# make nice output
	mu <- iTheta[elim.idx]
	Ynames <- rownames(sm_est$mu)
	dim(mu) <- c(length(mu),1)
	rownames(mu) <- Ynames
	dimnames(Ohat) <- list(Ynames,Ynames)
	retval <- list(mu=mu,Ohat=Ohat,n=sm_est$n,pp=pp)
	return(retval)
}

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
