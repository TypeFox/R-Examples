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
# Created: 2015.11.23
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

#' @include AllClass.r
#' @include utils.r
NULL

#' @title Matrix-wise Multivariate Operations
#'
#' @description 
#'
#' Element-wise multivariate operations. 
#'
#' @details
#'
#' These operations are operations on matrices: compute the symmetric square root
#' or the Cholesky factor. In the future, the matrix exponent and logarithm may
#' be implemented?
#'
#' @include AllClass.r
#' @inheritParams expm::sqrtm
#' @param x \code{madness} object.
#' @param ... further arguments passed to or from other methods.
#' @name matwise
#' @template etc
NULL

# 2FIX: add logm, expm

#' @name matwise
#' @rdname matwise
#' @aliases sqrtm
#' @exportMethod sqrtm
setGeneric('sqrtm', signature="x", function(x) standardGeneric('sqrtm'))
#' @rdname matwise
#' @aliases sqrtm,madness-method
setMethod("sqrtm", signature(x="madness"),
					function(x) {
						xtag <- x@xtag
						val <- expm::sqrtm(x@val)
						scalby <- (t(val) %x% diag(dim(val)[2])) + (diag(dim(val)[2]) %x% val)
						dvdx <- solve(scalby,x@dvdx)
						vtag <- paste0('sqrtm(',x@vtag,')')
						varx <- x@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})

#' @rdname matwise
#' @export 
#' @method chol madness
#' @usage chol(x,...)
#' @aliases chol
chol.madness <- function(x,...) {
 	xtag <- x@xtag
 	val <- chol(x@val)
	scalby <- t(val) %x% diag(dim(val)[2]) 
	scalby <- scalby + .do_commutator(t(val),scalby)
	Lm <- matrixcalc::L.matrix(dim(val)[2])
	scalby <- Lm %*% (scalby %*% t(Lm))
	dvdx <- .do_commutator(val,t(Lm) %*% solve(scalby,Lm %*% x@dvdx))
	vtag <- paste0('chol(',x@vtag,')')
	varx <- x@varx

	new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
}

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
