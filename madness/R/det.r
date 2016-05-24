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
# Created: 2015.11.18
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

#' @include AllClass.r
#' @include utils.r
NULL

#det <- function(x, ...) UseMethod('det')

# see also http://stackoverflow.com/a/12101238/164611
# for the mess of S3 S4 and S3.S4
#' @title Matrix Determinant
#'
#' @description
#'
#' Compute the determinant of a matrix. As for \code{base::determinant},
#' a list of the modulus and sign are returned.
#'
#'
#' @include AllClass.r
#' @param x \code{madness} object.
#' @inheritParams base::determinant
#' @note throws an error for non-square matrices or non-matrix input.
#' @return a list with elements \code{modulus} and \code{sign}, 
#' which are \code{madness} objects.
#' @name det
#' @template etc
NULL

#determinant <- function(x, logarithm=TRUE, ...) UseMethod('determinant')

# @method determinant madness
# @aliases determinant
# det#FOLDUP
# don't repeat yourself
#' @rdname det
#' @method determinant madness
#' @export
determinant.madness <- function(x,logarithm=TRUE,...) {
 	xtag <- x@xtag
	vtag <- paste0('determinant(',x@vtag,', logarithm=',as.character(logarithm),')')
	vdet <- determinant(x@val,logarithm=logarithm)

	varx <- x@varx

	# the value first
	val <- matrix(vdet$modulus)
	dvdx <- matrix(t(solve(x@val)),nrow=1) %*% x@dvdx
	if (!logarithm) { 
		dvdx <- val[1,1,drop=TRUE] * dvdx
	}

	retv <- list()
	retv$modulus <- new("madness", val=val, dvdx=dvdx, vtag=paste0(vtag,'$modulus'), xtag=xtag, varx=varx)

	# now the sign
	val <- matrix(vdet$sign)
	dvdx <- ifelse(val[1,1,drop=TRUE] == 0,NaN,0.0)
	dvdx <- matrix(dvdx,nrow=1,ncol=ncol(x@dvdx))

	retv$sign <- new("madness", val=val, dvdx=dvdx, vtag=paste0(vtag,'$sign'), xtag=xtag, varx=varx)

	invisible(retv)
}
# det <- function(x,...) { base::det(x,...) }
#' @rdname det
#' @export 
det <- base::det
environment(det) <- environment()

#' @rdname det
#' @exportMethod determinant
#' @aliases determinant
setGeneric('determinant', function(x,logarithm=TRUE,...) standardGeneric('determinant'))

#' @rdname det
#' @aliases determinant,madness,missing-method
setMethod("determinant", signature(x="madness"),determinant.madness)

#' @rdname det
#' @aliases determinant,madness,missing-method
setMethod("determinant", signature(x="madness",logarithm="missing"),determinant.madness)

#' @rdname det
#' @aliases determinant,madness,logical-method
setMethod("determinant", signature(x="madness",logarithm="logical"),determinant.madness)

#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
