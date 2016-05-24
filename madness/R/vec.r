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

#' @title vectorize a multidimensional array.
#'
#' @description 
#'
#' Turn a multidimensional array into a (column) vector. 
#' Turn a (typically symmetric) matrix into a (column) vector of
#' the lower triangular part. Or reverse these
#' operations.
#'
#' @param x a \code{madness} object or multidimensional array or matrix.
#' @param k the diagonal from which to subselect. 
#' @param symmetric logical whether to put the array on the antidiagonal
#' as well. Will throw an error if \code{k > 0}.
#' @return a \code{madness} object or an array, of the vectorized array
#' or the subselected part. For the inverse operations, promotes to a
#' \code{madness} of a matrix, or a matrix.
#' @seealso \code{\link{reshapes}}, in particular \code{tril}.
#' @template etc
#' @rdname vec
#' @name vec
#' @examples 
#' y <- matrix(rnorm(16),ncol=4)
#' sy <- y + t(y)
#' vy <- vec(sy)
#' vmy <- vec(madness(sy))
#' vhy <- vech(sy)
#' vmhy <- vech(madness(sy))
#'
#' ivech(c(1,2,3))
#' ivech(c(1,2,3),-1)
#' ivech(c(1,2,3),-1,symmetric=TRUE)
#' ivech(c(1,2,3,4,5,6,7,8),1)
NULL

#' @rdname vec
#' @aliases vec
#' @exportMethod vec
setGeneric('vec', function(x) standardGeneric('vec'))
#' @rdname vec
#' @aliases vec,madness-method
setMethod("vec", signature(x="madness"),
					function(x) {
						xtag <- x@xtag
						val <- x@val
						dim(val) <- c(prod(dim(val)),1)
						dvdx <- x@dvdx
						vtag <- paste0('vec(',x@vtag,')')
						varx <- x@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})
#' @rdname vec
#' @aliases vec,array-method
setMethod("vec", signature(x="array"),
					function(x) {
						dim(x) <- c(length(x),1)
						x
					})

#' @rdname vec
#' @aliases vech
#' @exportMethod vech
setGeneric('vech', function(x,k=0) standardGeneric('vech'))
#' @rdname vec
#' @aliases vech,array-method
setMethod("vech", signature(x="array"),
					function(x,k=0) {
						takeus <- row(x) >= col(x) - k
						x <- x[takeus]
						dim(x) <- c(length(x),1)
						x
					})

#' @rdname vec
#' @aliases vech,madness-method
setMethod("vech", signature(x="madness"),
					function(x,k=0) {
						xtag <- x@xtag
						val <- x@val
						takeus <- row(val) >= col(val) - k
						val <- val[takeus]
						dim(val) <- c(length(val),1)
						dvdx <- x@dvdx[which(takeus),]
						vtag <- paste0('vech(',x@vtag,')')
						varx <- x@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})

#' @rdname vec
#' @aliases ivech
#' @exportMethod ivech
setGeneric('ivech', function(x,k=0,symmetric=FALSE) standardGeneric('ivech'))

#' @rdname vec
#' @aliases ivech,ANY-method
setMethod("ivech", signature(x="ANY"),
					function(x,k=0,symmetric=FALSE) {
						len <- length(x)
						stopifnot((k <= 0) || (!symmetric))
						if (k <= 0) {
							Ms <- .quadeq(1,2*k+1,k*(k+1) - 2 * len)
						} else {
							Ms <- .quadeq(1,2*k+1,-k*(k+1) - 2 * len)
						}
						isok <- ((abs(Ms %% 1) <= 1e-12) & (Ms > 0))
						stopifnot(any(isok))
						M <- Ms[isok][1]
						retv <- array(0,dim=c(M,M))
						putus <- row(retv) >= col(retv) - k
						retv[putus] <- x
						if (symmetric) {
							# possibly double put on diagonal and not efficient, but no worries...
							retv <- t(retv)
							retv[putus] <- x
							retv <- t(retv)
						}
						retv
					})

#' @rdname vec
#' @aliases ivech,madness-method
setMethod("ivech", signature(x="madness"),
					function(x,k=0,symmetric=FALSE) {
						len <- length(x@val)
						stopifnot((k <= 0) || (!symmetric))
						if (k <= 0) {
							Ms <- .quadeq(1,2*k+1,k*(k+1) - 2 * len)
						} else {
							Ms <- .quadeq(1,2*k+1,-k*(k+1) - 2 * len)
						}
						isok <- ((abs(Ms %% 1) <= 1e-12) & (Ms > 0))
						stopifnot(any(isok))
						M <- Ms[isok][1]
						val <- array(0,dim=c(M,M))
						dvdx <- matrix(0,nrow=length(val),ncol=ncol(x@dvdx))
						putus <- row(val) >= col(val) - k
						val[putus] <- x@val
						dvdx[as.numeric(which(putus)),] <- x@dvdx
						if (symmetric) {
							# possibly double put on diagonal and not efficient, but no worries...
							val <- t(val)
							val[putus] <- x@val
							dvdx <- .do_commutator(val,dvdx)
							dvdx[as.numeric(which(putus)),] <- x@dvdx
							val <- t(val)
							dvdx <- .do_commutator(val,dvdx)
						}
						xtag <- x@xtag
						vtag <- paste0('ivech(',x@vtag,', ',k,', ',symmetric,')')
						varx <- x@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
