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
# Created: 2015.11.22
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

#' @include AllClass.r
#' @include utils.r
NULL

#' @title Sum and Product.
#'
#' @description
#'
#' Compute sum or product of \code{madness} objects.
#'
#' @include AllClass.r
#' @param x a numeric or \code{madness} object.
#' @param na.rm logical. Should missing values (including \sQuote{NaN}) be
#' removed?
#' @param ... ignored here.
#' @inheritParams base::sum
#' @return a \code{madness} object representing a scalar value.
#' @name sumprod
#' @template etc
#' @examples 
#' xv <- matrix(rnorm(5*5),ncol=5)
#' xmad <- madness(xv)
#' prod(xv)
#' sum(xv)
NULL
#' @rdname sumprod
#' @aliases sum
#' @aliases sum,madness-class
setMethod("sum", signature(x="madness"),
					function(x, ..., na.rm = FALSE) {
						xtag <- x@xtag
						val <- sum(x@val,na.rm=na.rm)

						if (na.rm) { 
							isok <- !(is.na(x@val) | is.nan(x@val))
							dvdx <- colSums(x@dvdx[which(isok),,drop=FALSE])
						} else {
							dvdx <- colSums(x@dvdx)
						}
						dvdx <- matrix(dvdx,nrow=1)

						vtag <- paste0('sum(',x@vtag,', na.rm=',na.rm,')')
						varx <- x@varx

						retv <- new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
						retv <- retv + sum(...,na.rm=na.rm)
						retv
					})
#' @rdname sumprod
#' @aliases prod
#' @aliases prod,madness-class
setMethod("prod", signature(x="madness"),
					function(x, ..., na.rm = FALSE) {
						xtag <- x@xtag
						if (na.rm) { 
							valiii <- which(!(is.na(x@val) | is.nan(x@val)))
						} else {
							valiii <- seq_len(length(x@val))
						}
						cval <- c(1,cumprod(x@val[valiii]))
						val <- cval[length(cval)]
						dvdx <- matrix(0,nrow=1,ncol=ncol(x@dvdx))
						for (jjj in seq_along(valiii)) {
							vi <- valiii[jjj]
							dvdx <- x@val[vi] * dvdx + cval[jjj] * x@dvdx[vi,,drop=FALSE]
						}
						dvdx <- matrix(dvdx,nrow=1)

						vtag <- paste0('prod(',x@vtag,', na.rm=',na.rm,')')
						varx <- x@varx

						retv <- new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
						retv <- retv * prod(...,na.rm=na.rm)
						retv
					})

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
