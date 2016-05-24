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
# Created: 2015.12.13
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav
# /usr/bin/r

#' @include AllClass.r
#' @include utils.r
NULL

#' @title Numerical (approximate) Differentiation.
#'
#' @description
#'
#' Approximates the derivative of a function at the input
#' by numerical methods. 
#'
#' @details
#'
#' For a multivariate-valued function of multivariate
#' data, approximates the derivative at a point via the
#' forward, central, or backward first differences, 
#' returning a \code{madness} object.
#'
#' @usage
#'
#' numderiv(f, x, eps=1e-8, type=c('forward','central','backward'),...)
#'
#' @param f a function, to be evaluated at and near \code{x}.
#' @param x array, matrix, or \code{madness} object.
#' @param eps the 'epsilon', a small value added or subtracted from \code{x} to
#' compute the first differences.
#' @param type the type of first difference, case-insensitive, substrings ok.
#' @param ... arguments passed on to \code{f}.
#' @return A matrix if \code{x} is numeric; a \code{madness} object if
#' \code{x} is a \code{madness} object.
#' @template etc
#' @examples 
#' f <- function(x,h) {
#'   cos(x + h)
#' }
#' x <- array(rnorm(100),dim=c(10,10))
#' madx <- numderiv(f,x,1e-8,h=pi)
#' @export
#' @exportMethod numderiv
#' @name numderiv
setGeneric('numderiv', function(f,x,eps=1e-8,type=c('forward','central','backward'),...) standardGeneric('numderiv'))
#' @rdname numderiv
#' @aliases ANY,array-method
setMethod('numderiv', signature(f='ANY', x='array'), 
	function(f,x,eps=1e-8,type=c('forward','central','backward'),...) {
	type <- match.arg(type)
	yval <- f(x,...)
	dapx <- matrix(0,length(yval),length(x))
	for (iii in seq_len(length(x))) {
		dydx <- switch(substring(tolower(type),1,2),
			fo={ 
				xalt <- x
				xalt[iii] <- xalt[iii] + eps
				yplus <- f(xalt,...)
				(yplus - yval) / eps },
			ba={ 
				xalt <- x
				xalt[iii] <- xalt[iii] - eps
				yneg <- f(xalt,...)
				(yval - yneg) / eps },
			ce={
				xalt <- x
				xalt[iii] <- xalt[iii] + eps
				yplus <- f(xalt,...)
				xalt <- x
				xalt[iii] <- xalt[iii] - eps
				yneg <- f(xalt,...)
				(yplus - yneg) / (2*eps) }
		)
		dapx[,iii] <- as.numeric(dydx)
	}
	dapx
})
#' @rdname numderiv
#' @aliases ANY,madness-method
setMethod('numderiv', signature(f='ANY', x='madness'), 
	function(f,x,eps=1e-8,type=c('forward','central','backward'),...) {
		type <- match.arg(type)
		val <- f(x@val,...)
		dfdx <- numderiv(f,x@val,eps=eps,type=type,...)
		dvdx <- dfdx %*% x@dvdx
		fname <- deparse(substitute(f))
		vtag <- paste0(fname,'(',x@vtag,')')
		xtag <- x@xtag
		varx <- x@varx
		new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
	})

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
