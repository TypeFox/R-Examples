# /usr/bin/r
#
# Copyright 2015-2016 Steven E. Pav. All Rights Reserved.
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
# Created: 2016.01.07
# Copyright: Steven E. Pav, 2015-2016
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

#' @include AllClass.r
#' @include utils.r
NULL

#' Diagonal Operations 
#'
#' @include AllClass.r
#' @param x \code{madness} object.
#' @param nrow,ncol Optional dimensions for the result when \code{x} is not a matrix value (NYI here)
#' @inheritParams base::diag
#' @note the (somewhat odd) use of \code{stats::diag} for two different
#' functions is \emph{not} repeated here, at least for now.
#' The \code{nrow} and \code{ncol} are ignored.
#' @seealso \code{\link{reshapes}}
#' @name todiag
#' @template etc
NULL

#' @rdname todiag
#' @aliases diag
#' @exportMethod diag
#' @inheritParams base::diag
setGeneric('diag', function(x,nrow,ncol) standardGeneric('diag'))
#' @rdname todiag
#' @aliases diag,madness-method
setMethod("diag", signature(x="madness"),
					function(x,nrow,ncol) {
						xtag <- x@xtag
						val <- x@val
						# 2FIX: this should just use the diag method, so that if e.g.
						# nrow and ncol are given, it does the right thing?
						takeus <- row(val) == col(val)
						val <- val[takeus]
						dim(val) <- c(length(val),1)
						dvdx <- x@dvdx[which(takeus),]
						vtag <- paste0('diag(',x@vtag,')')
						varx <- x@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})

#' @rdname todiag
#' @aliases todiag
#' @exportMethod todiag
setGeneric('todiag', function(x) standardGeneric('todiag'))
#' @rdname todiag
#' @aliases todiag,madness-method
setMethod("todiag", signature(x="madness"),
					function(x) {
						xtag <- x@xtag
						val <- diag(as.numeric(x@val))
						takeus <- row(val) == col(val)
						dvdx <- matrix(0,nrow=length(takeus),ncol=ncol(x@dvdx))
						dvdx[which(takeus),] <- x@dvdx
						vtag <- paste0('todiag(',x@vtag,')')
						varx <- x@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
