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

#' @title Element-wise Multivariate Operations
#'
#' @description 
#'
#' Element-wise multivariate operations. 
#'
#' @details
#'
#' These operations are scalar-to-scalar operations applied to each element of a 
#' multidimensional array. 
#'
#' @note
#'
#' The \code{exp}, \code{log}, and \code{sqrt}
#' functions are not to be confused with the matrix-wise operations,
#' \code{expm}, \code{logm} and \code{sqrtm}
#'
#' @include AllClass.r
#' @param x \code{madness} object.
#' @seealso \link{matwise}
#' @name elwise
#' @template etc
NULL

#' @rdname elwise
#' @aliases abs,madness-method
setMethod("abs", signature(x="madness"),
					function(x) {
						xtag <- x@xtag
						val <- abs(x@val)
						scald <- sign(x@val)
						scald[scald==0] <- NaN
						dvdx <- as.numeric(scald) * x@dvdx
						vtag <- paste0('abs(',x@vtag,')')
						varx <- x@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})

#' @rdname elwise
#' @aliases exp,madness-method
setMethod("exp", signature(x="madness"),
					function(x) {
						xtag <- x@xtag
						val <- exp(x@val)
						scald <- val
						dvdx <- as.numeric(scald) * x@dvdx
						vtag <- paste0('exp(',x@vtag,')')
						varx <- x@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})

#' @rdname elwise
#' @aliases log,madness-method
setMethod("log", signature(x="madness"),
					function(x) {
						xtag <- x@xtag
						val <- log(x@val)
						scald <- 1 / x@val
						dvdx <- as.numeric(scald) * x@dvdx
						vtag <- paste0('log(',x@vtag,')')
						varx <- x@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})

#' @rdname elwise
#' @aliases log10,madness-method
setMethod("log10", signature(x="madness"),
					function(x) {
						xtag <- x@xtag
						val <- log10(x@val)
						scald <- 1 / (log(10) * x@val)
						dvdx <- as.numeric(scald) * x@dvdx
						vtag <- paste0('log10(',x@vtag,')')
						varx <- x@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})

#' @rdname elwise
#' @aliases sqrt,madness-method
setMethod("sqrt", signature(x="madness"),
					function(x) {
						xtag <- x@xtag
						val <- sqrt(x@val)
						scald <- 1 / (2*val)
						dvdx <- as.numeric(scald) * x@dvdx
						vtag <- paste0('sqrt(',x@vtag,')')
						varx <- x@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})

#' @rdname elwise
#' @aliases sin,madness-method
setMethod("sin", signature(x="madness"),
					function(x) {
						xtag <- x@xtag
						val <- sin(x@val)
						scald <- cos(x@val)
						dvdx <- as.numeric(scald) * x@dvdx
						vtag <- paste0('sin(',x@vtag,')')
						varx <- x@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})

#' @rdname elwise
#' @aliases cos,madness-method
setMethod("cos", signature(x="madness"),
					function(x) {
						xtag <- x@xtag
						val <- cos(x@val)
						scald <- -sin(x@val)
						dvdx <- as.numeric(scald) * x@dvdx
						vtag <- paste0('cos(',x@vtag,')')
						varx <- x@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})

#' @rdname elwise
#' @aliases tan,madness-method
setMethod("tan", signature(x="madness"),
					function(x) {
						xtag <- x@xtag
						val <- tan(x@val)
						scald <- 1 + (val^2)
						dvdx <- as.numeric(scald) * x@dvdx
						vtag <- paste0('tan(',x@vtag,')')
						varx <- x@varx

						new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
					})




#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
