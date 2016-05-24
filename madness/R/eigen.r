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
# Created: 2016.01.09
# Copyright: Steven E. Pav, 2016
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

#' @include AllClass.r
#' @include utils.r
NULL

#' @title Spectral Decomposition of a Matrix
#'
#' @description 
#'
#' Computes eigenvalues and eigenvectors of numeric (double, integer, logical) or 
#' complex \code{madness} matrices.
#'
#' @details
#'
#' The singular value decomposition of the matrix \eqn{X}{X} is
#' \deqn{X = U D V',}{X = U D V',}
#' where \eqn{U} and \eqn{V} are orthogonal, \eqn{V'} is \eqn{V}
#' transposed, and \eqn{D} is a diagonal matrix with the singular
#' values on the diagonal.
#'
#' @include AllClass.r
#' @inheritParams base::eigen
#' @param x \code{madness} object representing a numeric matrix
#' whose spectral decomposition is to be computed.
#' @seealso \code{\link{eigen}}.
#' @return a list with components
#' \describe{
#' \item{values}{a \code{madness} object of a vector containing
#' the \eqn{p} eigenvalues of \code{x}, sorted in \emph{decreasing} order,
#' according to \code{Mod(value)} in the assymetric case when they might
#' be complex (even for real matrices). For real asymmetric matrices
#' the vector will be complex only if complex conjugate pairs of eigenvalues are 
#' detected.}
#' \item{vectors}{either a \eqn{p \times p}{p * p} matrix whose columns contain the
#' eigenvectors of \code{x} or \code{NULL} if \code{only.values} is 
#' \code{TRUE}. The vectors are normalized to unit length.
#'
#' Recall that the eigenvectors are only defined up to a constant: 
#' even when the length is specified they are still only defined up to a 
#' scalar of modulus one (the sign for real matrices).  
#' If \code{r <- eigen(A)}, and \code{V <- r$vectors; lam <- r$values}, then
#' \deqn{A = V Lmbd V^{-1}}{A = V Lmbd V^(-1)}
#' (up to numerical fuzz), where \code{Lmbd =diag(lam)}.
#' }
#' }
#'
#' @references
#'
#' Izenman, Alan Julian. "Reduced-Rank Regression for the Multivariate Linear
#' Model." Journal of Multivariate Analysis 5, pp 248-264 (1975).
#' \url{http://www.sciencedirect.com/science/article/pii/0047259X75900421}
#'
#' Kato, Tosio. "Perturbation Theory for Linear Operators."
#' Springer (1995).
#' \url{http://www.maths.ed.ac.uk/~aar/papers/kato1.pdf}
#'
#' @name eigen
#' @template etc
NULL

#' @name eigen
#' @rdname eigen
#' @exportMethod eigen
setGeneric('eigen', function(x,symmetric,only.values=FALSE,EISPACK=FALSE) standardGeneric('eigen'))
#' @rdname eigen
#' @aliases eigen,madness-method
setMethod("eigen", signature(x="madness"),
					function(x,symmetric,only.values=FALSE,EISPACK=FALSE) {
						if (!missing(symmetric) && symmetric) {
							# if symmetric is set true, then we use only the lower
							# triangular part? in which case we need to enforce symmetry...
							vtag <- x@vtag
							x <- vech(x,k=0)
							x <- ivech(x,k=0,symmetric=TRUE)
							x@vtag <- vtag
						}

						xtag <- x@xtag
						if (missing(symmetric)) {
							vvals <- eigen(x@val,only.values=only.values,EISPACK=EISPACK)
						} else {
							vvals <- eigen(x@val,symmetric=symmetric,only.values=only.values,EISPACK=EISPACK)
						}

						# eigen*values*
						p <- length(vvals$values)
						rowr <- lapply(seq_len(p),function(jj) { (t(vvals$vectors[,jj,drop=FALSE]) %x% t(vvals$vectors[,jj,drop=FALSE])) %*% x@dvdx })
						dvdx <- do.call(rbind,rowr)
						vtag <- paste0('eigen(',x@vtag,')$values')
						varx <- x@varx
						retv <- list()
						retv$values <- new("madness", val=array(vvals$values,dim=c(1,p)), dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
						
						# eigen*vectors* are trickier
						rowr <- lapply(seq_len(p),
													 function(jj) {
														 scals <- 1 / (vvals$values[jj] - vvals$values[-jj])
														 krons <- t(vvals$vectors[,-jj,drop=FALSE]) %x% t(vvals$vectors[,jj,drop=FALSE]) %*% x@dvdx
														 rv <- vvals$vectors[,-jj,drop=FALSE] %*% diag(scals) %*% krons
													 })
						dvdx <- do.call(rbind,rowr)
						vtag <- paste0('eigen(',x@vtag,')$vectors')
						varx <- x@varx
						retv$vectors <- new("madness", val=vvals$vectors, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
						retv
					})

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
