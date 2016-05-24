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

#' Matrix Trace
#'
#' @include AllClass.r
#' @param x \code{madness} object.
#' @exportMethod matrix.trace
#' @name matrix.trace
#' @template etc
setGeneric('matrix.trace', function(x) standardGeneric('matrix.trace'))

# det#FOLDUP

##' @rdname matrix.trace
##' @aliases matrix.trace,matrix-method
#setMethod("matrix.trace", signature(x="matrix"),function(x) { sum(diag(x)) })

#' @rdname matrix.trace
#' @aliases matrix.trace,ANY-method
setMethod("matrix.trace", signature(x="ANY"),
					function(x) { sum(diag(x)) })

#' @rdname matrix.trace
#' @aliases matrix.trace,madness-method
setMethod("matrix.trace", signature(x="madness"),
					function(x) {
						xtag <- x@xtag
						val <- x@val
						takeus <- row(val) == col(val)
						val <- array(sum(val[takeus]),dim=c(1,1))

						dvdx <- matrix(colSums(x@dvdx[which(takeus),,drop=FALSE]),nrow=1)
						vtag <- paste0('matrix.trace(',x@vtag,')')
						varx <- x@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})

#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
