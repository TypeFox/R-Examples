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

#' Basic Matrix Inversion 
#'
#' @include AllClass.r
#' @param a,b \code{madness} object or matrix value.
#' @inheritParams base::solve
#' @name solve
#' @template etc
#' @exportMethod solve
setGeneric('solve', function(a,b) standardGeneric('solve'))

# solve#FOLDUP

#' @rdname solve
#' @aliases solve,ANY,missing-method
setMethod("solve", signature(a="ANY",b="missing"),
					function(a,b) { base::solve(a) })

#' @rdname solve
#' @aliases solve,madness,missing-method
setMethod("solve", signature(a="madness",b="missing"),
					function(a,b) {
						xtag <- a@xtag
						val <- solve(a@val)

						dvdx <- - (t(val) %x% val) %*% a@dvdx
						vtag <- paste0('solve(',a@vtag,')')
						varx <- a@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})

#' @rdname solve
#' @aliases solve,madness,madness-method
setMethod("solve", signature(a="madness",b="madness"),
					function(a,b) {
						xtag <- .check_common_xtag(a,b)
						val <- solve(a@val,b@val)

						# wasteful in that we invert twice, but whatever.
						ainv <- solve(a@val)
						dvdx <- (diag(dim(val)[2]) %x% ainv) %*% b@dvdx -
							(t(val) %x% ainv) %*% a@dvdx
						vtag <- paste0('solve(',a@vtag,', ',b@vtag,')')
						varx <- .get_a_varx(a,b)

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})

#' @rdname solve
#' @aliases solve,madness,array-method
setMethod("solve", signature(a="madness",b="array"),
					function(a,b) {
						xtag <- a@xtag
						val <- solve(a@val,b)

						# wasteful in that we invert twice, but whatever.
						dvdx <- - (t(val) %x% solve(a@val)) %*% a@dvdx
						vtag <- paste0('solve(',a@vtag,', numeric)')
						varx <- a@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})

#' @rdname solve
#' @aliases solve,madness,ANY-method
setMethod("solve", signature(a="madness",b="ANY"),
					function(a,b) {
						xtag <- a@xtag
						val <- solve(a@val,b)

						# wasteful in that we invert twice, but whatever.
						dvdx <- - (t(val) %x% solve(a@val)) %*% a@dvdx
						vtag <- paste0('solve(',a@vtag,', numeric)')
						varx <- a@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})

#' @rdname solve
#' @aliases solve,array,madness-method
setMethod("solve", signature(a="array",b="madness"),
					function(a,b) {
						xtag <- b@xtag
						val <- solve(a,b@val)

						# wasteful in that we invert twice, but whatever.
						dvdx <- (diag(dim(val)[2]) %x% solve(a)) %*% b@dvdx 
						vtag <- paste0('solve(numeric, ',b@vtag,')')
						varx <- b@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})

#' @rdname solve
#' @aliases solve,ANY,madness-method
setMethod("solve", signature(a="ANY",b="madness"),
					function(a,b) {
						xtag <- b@xtag
						val <- solve(a,b@val)

						# wasteful in that we invert twice, but whatever.
						dvdx <- (diag(dim(val)[2]) %x% solve(a)) %*% b@dvdx 
						vtag <- paste0('solve(numeric, ',b@vtag,')')
						varx <- b@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})

#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
