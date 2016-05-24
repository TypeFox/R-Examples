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
# Created: 2015.11.23
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

#' @include AllClass.r
#' @include utils.r
NULL

#' @title Matrix and vector norms.
#'
#' @description 
#'
#' Compute the norm of a vector or matrix, as determined by the
#' \code{type}.
#'
#' @include AllClass.r
#' @inheritParams base::norm
#' @param x \code{madness} object.
#' @param ... further arguments passed to or from other methods.
#' @return the matrix norm, a non-negative number.
#' @note This should probably be fixed to return a scalar, not a 1 by 1 matrix?
#' @name norm
#' @template etc
NULL

#' @name norm
#' @rdname norm
#' @aliases maxeig
#' @exportMethod maxeig
setGeneric('maxeig', signature="x", function(x) standardGeneric('maxeig'))
#' @rdname norm
#' @aliases maxeig,madness-method
setMethod("maxeig", signature(x="madness"),
					function(x) {
						xtag <- x@xtag
						vvals <- svd(x@val,nu=1,nv=1)
						val <- matrix(vvals$d[1])
						dvdx <- (t(vvals$v) %x% t(vvals$u)) %*% x@dvdx
						vtag <- paste0('maxeig(',x@vtag,')')
						varx <- x@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})

.normit <- function(x,type='One') {
	type <- tolower(substr(type,0,1))
	stopifnot(type %in% c('o','1','2','i','f','m'))
	xtag <- x@xtag
	#val <- norm(x@val,type=type)
	nr <- nrow(x@val)
	nc <- ncol(x@val)
	if (type %in% c('o','1')) {
# DRY: colsums
		cs <- colSums(abs(x@val))
		midx <- which.max(cs)
		val <- cs[midx]
		dvdx <- array(as.numeric(sign(x@val[,midx])),dim=c(1,nr)) %*% x@dvdx[which(col(x@val)==midx),,drop=FALSE]
	} else if (type == 'i') {
# DRY: rowsums
		rs <- rowSums(abs(x@val))
		midx <- which.max(rs)
		val <- rs[midx]
		dvdx <- array(as.numeric(sign(x@val[midx,])),dim=c(1,nc)) %*% x@dvdx[which(row(x@val)==midx),,drop=FALSE]
	} else if (type == 'f') {
# DRY:
		val <- sqrt(sum(x@val^2))
		dvdx <- (1/val) * (array(x@val,dim=c(1,length(x@val))) %*% x@dvdx)
	} else if (type == 'm') {
		val <- max(abs(x@val))
		midx <- which.max(abs(x@val))
		dvdx <- sign(x@val[midx]) * x@dvdx[midx,,drop=FALSE]
	} else {
		# DRY: farm this off to maxeig?
		vvals <- svd(x@val,nu=1,nv=1)
		val <- matrix(vvals$d[1])
		dvdx <- (t(vvals$v) %x% t(vvals$u)) %*% x@dvdx
	}

	vtag <- paste0('norm(',x@vtag,", '",type,"')")
	varx <- x@varx
	val <- array(val,dim=c(1,1))
	dvdx <- array(dvdx,dim=c(1,length(dvdx)))

	new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
}

#' @name norm
#' @rdname norm
#' @aliases norm
#' @exportMethod norm
setGeneric('norm', function(x,type) standardGeneric('norm'))
#' @rdname norm
#' @aliases maxeig,madness-method
setMethod("norm", signature(x="madness",type='missing'), function(x) .normit(x))
#' @rdname norm
#' @aliases maxeig,madness-method
setMethod("norm", signature(x="madness",type='ANY'), .normit)

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
