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

#' Basic Reshape Operations
#'
#' @include AllClass.r
#' @param x \code{madness} object.
#' @param value an array of the new dimensions of the object value.
#' @inheritParams Matrix::tril
#' @param k the index of the diagonal number from which to extract.\code{madness} object.
#' @seealso \code{\link{vec}}, \code{\link{todiag}}
#' @name reshapes
#' @template etc
NULL

# transpose#FOLDUP

#' @rdname reshapes
#' @aliases t
#' @exportMethod t
setGeneric('t', function(x) standardGeneric('t'))
#' @rdname reshapes
#' @aliases t,madness-method
setMethod("t", signature(x="madness"),
					function(x) {
						xtag <- x@xtag
						val <- t(x@val)
						dvdx <- .do_commutator(val,x@dvdx)
						vtag <- paste0('t(',x@vtag,')')
						varx <- x@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})

#' @rdname reshapes
#' @aliases tril
#' @exportMethod tril
setGeneric('tril', function(x,k=0,...) standardGeneric('tril'))
# 2FIX: must I check if this has already been defined as a generic?
#' @rdname reshapes
#' @aliases tril,madness-method
setMethod("tril", signature(x="madness"),
					function(x,k=0) {
						xtag <- x@xtag
						val <- x@val
						takeus <- row(val) >= col(val) - k
						val[!takeus] <- 0
						dvdx <- x@dvdx
						dvdx[which(!takeus),] <- 0
						vtag <- paste0('tril(',x@vtag,', ',k,')')
						varx <- x@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})

#' @rdname reshapes
#' @aliases triu
#' @exportMethod triu
setGeneric('triu', function(x,k=0,...) standardGeneric('triu'))
#' @rdname reshapes
#' @aliases triu,madness-method
setMethod("triu", signature(x="madness"),
					function(x,k=0) {
						xtag <- x@xtag
						val <- x@val
						takeus <- row(val) <= col(val) - k
						val[!takeus] <- 0
						dvdx <- x@dvdx
						dvdx[which(!takeus),] <- 0
						vtag <- paste0('triu(',x@vtag,', ',k,')')
						varx <- x@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})


#' @rdname reshapes
#' @aliases dim<-,madness,ANY-method
setMethod("dim<-", signature(x="madness",value="ANY"),
					function(x,value) {
						xtag <- x@xtag
						val <- x@val
						dim(val) <- value
						dvdx <- x@dvdx
						vtag <- paste0('reshape(',x@vtag,', ',
													 as.character(enquote(value))[2],
													 ')')
						varx <- x@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})

# see http://stackoverflow.com/a/8057007/164611
#' Extract parts of a \code{madness} value.
#'
#' @param x a \code{madness} object.
#' @param j,...  further indices specifying elements to extract or 
#' replace.  
#' @inheritParams base::`[`
#' @name [
#' @aliases [,madness,ANY,ANY,ANY-method
#' @docType methods
#' @rdname extract-methods
#' @template etc
setMethod("[", signature(x="madness",i="ANY",j="ANY"),
					function(x,i,j,...,drop) {
						getidx <- vector(length=length(x@val))
						dim(getidx) <- dim(x@val)
						getidx[i,j,...] <- TRUE
						val <- x@val[i,j,...,drop=FALSE]
						dvdx <- x@dvdx[which(getidx),,drop=FALSE]

						retv <- new("madness", val=val, dvdx=dvdx,
												xtag=x@xtag,
												vtag=paste0(x@vtag,'[...]'),
												varx=x@varx)
					})

# 2FIX: define the setter method? but, wait we don't want to be
# able to poke elements willy nilly, right?

#UNFOLD

#' @rdname reshapes
#' @export 
#' @method aperm madness
#' @usage aperm(a, perm=NULL, resize=TRUE, ...)
#' @inheritParams base::aperm
#' @aliases aperm
aperm.madness <- function(a, perm=NULL, resize=TRUE, ...) {
 	xtag <- a@xtag
 	val <- aperm(a@val,perm=perm,resize=resize)
	oldids <- array(1:length(a@val),dim=dim(a@val))
	prmids <- aperm(oldids,perm=perm,resize=resize)
	dvdx <- a@dvdx[as.numeric(prmids),,drop=FALSE]
	vtag <- paste0('aperm(',a@vtag,', ',perm,')')
	varx <- a@varx

	new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
}

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
