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
# Created: 2015.11.20
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

#' @include AllClass.r
#' @include utils.r
NULL

#' Row and Column Bind
#'
#' @template etc
#' @include AllClass.r
#' @param x,y \code{madness} or array, numeric, matrix objects.
#' @param ... optional arguments for methods (ignored here).
#' @inheritParams base::cbind2
#' @name bind
NULL

# bind#FOLDUP

#' @rdname bind
#' @method c madness
#' @export
#' @usage \\method{c}{madness}(...)
c.madness <- function(...) { rbind(...) }  # nocov

# c.f. http://stackoverflow.com/a/28126631/164611

#' @rdname bind
#' @aliases cbind2,madness,missing-method
#' @exportMethod cbind2
setMethod("cbind2", signature(x="madness",y="missing"),
					function(x,y,...) {
						xtag <- x@xtag
						val <- x@val
						dvdx <- x@dvdx

						vtag <- paste0('cbind(',x@vtag,')')
						varx <- x@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})
#' @rdname bind
#' @aliases cbind2,madness,madness-method
setMethod("cbind2", signature(x="madness",y="madness"),
					function(x,y,...) {
						xtag <- .check_common_xtag(x,y)
						val <- cbind(x@val,y@val)
						dvdx <- rbind(x@dvdx,y@dvdx)

						vtag <- paste0('cbind(',x@vtag,',',y@vtag,')')
						varx <- .get_a_varx(x,y)

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})

#' @rdname bind
#' @aliases cbind2,madness,ANY-method
setMethod("cbind2", signature(x="madness",y="ANY"),
					function(x,y,...) {
						xtag <- x@xtag
						val <- cbind(x@val,y)
						dvdx <- rbind(x@dvdx,
													array(0,dim=c(length(y),ncol(x@dvdx))))

						vtag <- paste0('cbind(',x@vtag,', numeric)')
						varx <- x@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})
#' @rdname bind
#' @aliases cbind2,ANY,madness-method
setMethod("cbind2", signature(x="ANY",y="madness"),
					function(x,y,...) {
						xtag <- y@xtag
						val <- cbind(x,y@val)
						dvdx <- rbind(array(0,dim=c(length(x),ncol(y@dvdx))),
													y@dvdx)

						vtag <- paste0('cbind(numeric,',y@vtag,')')
						varx <- y@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})

# this is way trickier, unfortunately, because of numerator layout.
# so we apply the commutator to the dvdx, then rbind them, then apply 
# commutator again

#' @rdname bind
#' @aliases rbind2,madness,missing-method
#' @exportMethod rbind2
setMethod("rbind2", signature(x="madness",y="missing"),
					function(x,y,...) {
						xtag <- x@xtag
						val <- x@val
						dvdx <- x@dvdx

						vtag <- paste0('rbind(',x@vtag,')')
						varx <- x@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})
#' @rdname bind
#' @aliases rbind2,madness,madness-method
setMethod("rbind2", signature(x="madness",y="madness"),
					function(x,y,...) {
						xtag <- .check_common_xtag(x,y)
						val <- rbind(x@val,y@val)

						dvdx <- rbind(.do_commutator(t(x@val),x@dvdx),
													.do_commutator(t(y@val),y@dvdx))
						dvdx <- .do_commutator(val,dvdx)

						vtag <- paste0('rbind(',x@vtag,',',y@vtag,')')
						varx <- .get_a_varx(x,y)

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})

#' @rdname bind
#' @aliases rbind2,madness,ANY-method
setMethod("rbind2", signature(x="madness",y="ANY"),
					function(x,y,...) {
						xtag <- x@xtag
						val <- rbind(x@val,y)

						dvdx <- rbind(.do_commutator(t(x@val),x@dvdx),
													array(0,dim=c(length(y),ncol(x@dvdx))))
						dvdx <- .do_commutator(val,dvdx)

						vtag <- paste0('rbind(',x@vtag,',numeric)')
						varx <- x@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})

#' @rdname bind
#' @aliases rbind2,ANY,madness-method
setMethod("rbind2", signature(x="ANY",y="madness"),
					function(x,y,...) {
						xtag <- y@xtag
						val <- rbind(x,y@val)

						dvdx <- rbind(array(0,dim=c(length(x),ncol(y@dvdx))),
													.do_commutator(t(y@val),y@dvdx))
						dvdx <- .do_commutator(val,dvdx)

						vtag <- paste0('rbind(numeric,',y@vtag,')')
						varx <- y@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})

#UNFOLD

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
