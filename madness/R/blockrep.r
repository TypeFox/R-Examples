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
# Created: 2015.12.05
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

#' @include AllClass.r
#' @include utils.r
NULL

#' @title Replicate blocks of multidimensional value.
#'
#' @description 
#'
#' Replicates a multidimensional object a number of times along
#' given dimensions.
#'
#' @details
#'
#' Given a k-dimensional object, and an l-vector of positive
#' integers, for l >= k, copy the input object l_i times in
#' the ith dimension. Useful for replication and (slow, fake)
#' outer products.
#'
#' \code{repto} replicates to the given dimension, assuming the
#' given dimension are integer multiples of the input dimensions.
#'
#' @usage
#'
#' blockrep(x, nreps)
#'
#' repto(x, newdim)
#'
#' @param x a \code{madness} object, representing a k-dimensional object.
#' @param nreps an l-vector of positive integers, representing how
#' many times to copy the object.
#' @param newdim an l-vector of positive integers of the new dimension
#' of the output object. These must be integer multiples of the 
#' input dimensions.
#' @return A \code{madness} object replicated out.
#' @note
#' An error will be thrown if \code{nreps} or \code{newdim} are improper.
#' @template etc
#' @examples 
#' set.seed(123)
#' y <- array(rnorm(3*3),dim=c(3,3))
#' dy <- matrix(rnorm(length(y)*2),ncol=2)
#' dx <- crossprod(matrix(rnorm(ncol(dy)*100),nrow=100))
#' obj0 <- madness(val=y,vtag='y',xtag='x',dvdx=dy,varx=dx)
#'
#' anobj <- blockrep(obj0,c(1,2,1))
#' anobj <- blockrep(obj0,c(1,1,2))
#' anobj <- repto(obj0,c(9,12,4))
#' @export
#' @rdname blockrep
blockrep <- function(x, nreps) {
	olddim <- dim(x@val)
	if (length(nreps) < length(olddim)) {
		nreps <- c(nreps,rep(1,length(olddim)-length(nreps)))
	}
	stopifnot(length(nreps) >= length(olddim),
						all(nreps >= 1),
						all(abs(nreps %% 1) < 1e-12))
	nreps <- as.integer(nreps)
	xval <- x@val
	if (length(nreps) > length(olddim)) { 
		# so. stupid.
		dim(xval) <- c(dim(xval),rep(1,length(nreps) - length(olddim)))
	}
	rm(olddim)
	outdm <- length(nreps)

	# now replicate!
	# erk? http://stackoverflow.com/a/13034947/164611
	totd <- seq_len(outdm)
	xidx <- array(seq_len(length(xval)),dim=dim(xval))
	for (iii in which(nreps > 1)) {
		marg <- setdiff(totd,iii)
		prmd <- rep(1,outdm)
		prmd[c(iii,marg)] <- seq_len(outdm)

		xval <- aperm(apply(xval,MARGIN=marg,FUN=rep,nreps[iii]),prmd)
		xidx <- aperm(apply(xidx,MARGIN=marg,FUN=rep,nreps[iii]),prmd)
	}
	dvdx <- x@dvdx[as.numeric(xidx),,drop=FALSE]
	vtag <- paste0('blockrep(',x@vtag,', ', as.character(enquote(nreps))[2],')')
	xtag <- x@xtag
	varx <- x@varx

	new("madness", val=xval, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
}
#' @rdname blockrep
#' @export
repto <- function(x, newdim) {
	olddim <- dim(x@val)
	stopifnot(length(newdim) >= length(olddim),
						all(newdim >= 1))
	olddim <- c(olddim,rep(1,length(newdim)-length(olddim)))
	nreps <- newdim / olddim
	retv <- blockrep(x,nreps)
	retv@vtag <- paste0('repto(',x@vtag,', ', as.character(enquote(newdim))[2],')')
	retv
}

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
