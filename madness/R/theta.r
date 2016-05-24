# /usr/bin/r
#
# Copyright 2015-2015 Steven E. Pav. All Rights Reserved.
# Author: Steven E. Pav
#
# This file is part of madness.
#
# madness is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# madness is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with madness.  If not, see <http://www.gnu.org/licenses/>.
#
# Created: 2015.12.10
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

#' @include AllClass.r
#' @include utils.r
NULL

#' @title Estimate the symmetric second moment array of values.
#'
#' @description 
#'
#' Given rows of observations of some vector (or multidimensional
#' data), estimates the second moment by taking a simple mean,
#' returning a \code{madness} object.
#'
#' @details
#'
#' Given a 
#' \eqn{n\times k_1 \times k_2 \times ... \times k_l}{n x k_1 x k_2 ... x k_l}
#' array whose 'rows' are independent observations of \eqn{X}{X}, computes the 
#' \eqn{k_1 \times k_2 \times ... \times k_l \times k_1 \times k_2 ... k_l}{k_1 x k_2 x ... x k_l x k_1 x k_2 ... k_l}
#' array of the mean of \eqn{\mathrm{outer}(X,X)}{outer(X,X)} based on \eqn{n}{n} observations,
#' returned as a \code{\link{madness}} object. The variance-covariance
#' is also estimated, and stored in the object.
#'
#' One may use the default method for computing covariance,
#' via the \code{\link{vcov}} function, or via a 'fancy' estimator,
#' like \code{sandwich:vcovHAC}, \code{sandwich:vcovHC}, \emph{etc.}
#'
#' @usage
#'
#' theta(X, vcov.func=vcov, xtag=NULL)
#'
#' @param X a multidimensional array (or a data frame) of observed values.
#' @param vcov.func a function which takes an object of class \code{lm},
#' and computes a variance-covariance matrix. If equal to the string
#' "normal", we assume multivariate normal returns.
#' @param xtag an optional string tag giving the name of the input data.
#' defaults to figuring it out from the input expression.
#' @return A \code{madness} object representing the mean of the outer
#' product of the tail dimensions of \code{X}.
#' @template etc
#' @seealso \code{\link{twomoments}}
#' @examples 
#' set.seed(123)
#' X <- matrix(rnorm(1000*3),ncol=3)
#' th <- theta(X)
#'
#' \dontrun{
#' if (require(sandwich)) {
#'  th2 <- theta(X,vcov.func=vcovHC)
#' }
#' }
#' # works on data frames too:
#' set.seed(456)
#' X <- data.frame(a=runif(100),b=rnorm(100),c=1)
#' th <- theta(X)
#' @export
#' @rdname theta
theta <- function(X,vcov.func=vcov,xtag=NULL) {
	if (missing(xtag)) {
		xtag <- deparse(substitute(X))
	}
	X <- na.omit(X)
	if (is.data.frame(X)) {
		X <- as.matrix(X)
	}
	dimX <- dim(X)
	n <- dimX[1]
	p <- prod(dimX[2:length(dimX)])
	dim(X) <- c(n,p)
	# prealloc
	Y <- matrix(NA,nrow=n,ncol=p*(p+1)/2)
	offs <- 0
	for (iii in 1:p) {
		Y[,offs+1+(0:(p-iii))] = X[,iii] * X[,iii:p]
		offs <- offs + (p-iii) + 1
	}
	mod2 <- lm(Y ~ 1)
	rm(Y)

	mu <- mod2$coefficients
	dim(mu) <- c(length(mu),1)
	Ohat = vcov.func(mod2)
	rm(mod2)

	retv <- madness(val=mu, vtag='theta', xtag=xtag, varx=Ohat)
	# this has been vech'ed, so unvech it:
	retv <- ivech(retv,0,symmetric=TRUE)
	if (length(dimX) <= 2) {
		dim(retv) <- c(dimX[2],dimX[2])
	} else {
		dim(retv) <- c(dimX[2:length(dimX)],dimX[2:length(dimX)])
	}
	retv@vtag <- 'theta'
	retv
}

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
