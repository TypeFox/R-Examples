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
# Created: 2016.01.12
# Copyright: Steven E. Pav, 2016
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

#' @include AllClass.r
#' @include utils.r
NULL

#' @title Maxima and Minima
#'
#' @description 
#'
#' Return the maxima and minima of the input values.
#'
#' @details
#'
#' \code{max} and \code{min} return the maximum or minimum of \emph{all} the
#' values present in their arguments.
#'
#' If \code{na.rm} is \code{FALSE} and \code{NA} value in any of the arguments
#' will cause a value of \code{NA} to be returned, otherwise \code{NA} values are
#' ignored.
#'
#' The minimum and maximum of a numeric empty set are \code{+Inf} and 
#' \code{-Inf} (in this order!) which ensures \emph{transitivity}, e.g.,
#' \code{min(x1, min(x2)) == min(x1, x2)}.  For numeric \code{x} 
#' \code{max(x) == -Inf} and \code{min(x) == +Inf} whenever 
#' \code{length(x) == 0} (after removing missing values if requested).  
#'
#' @include AllClass.r
#' @inheritParams base::max
#' @param x \code{madness} object arguments.
#' @param ... \code{madness} object arguments.
#' @name max
#' @template etc
NULL

#' @rdname max
#' @aliases max
#' @aliases max,madness-class
setMethod("max", signature(x="madness"),
					function(x, ..., na.rm = FALSE) {
						xtag <- x@xtag
						val <- max(x@val,na.rm=na.rm)

						if (na.rm) { 
							dvdx <- x@dvdx[which.max(x@val),,drop=FALSE]
						} else {
							isok <- !(is.na(x@val) | is.nan(x@val))
							if (all(isok)) {
								dvdx <- x@dvdx[which.max(x@val),,drop=FALSE]
							} else {
								dvdx <- matrix(NA,nrow=1,ncol=ncol(x@dvdx))
							}
						}
						dvdx <- matrix(dvdx,nrow=1)

						vtag <- paste0('max(',x@vtag,', na.rm=',na.rm,')')
						varx <- x@varx

						retv <- new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
						if (length(list(...)) > 0) {
							altv <- max(...,na.rm=na.rm)

							# combine them:
							if (is(altv,'madness')) {
								xtag <- .check_common_xtag(retv,altv)
								varx <- .get_a_varx(retv,altv)
								vtag <- paste0('max(',retv@vtag,', ',altv@vtag,', na.rm=',na.rm,')')
								if (na.rm) {
									if (is.na(retv@val) || is.nan(retv@val) || (!is.na(altv@val) && !is.nan(altv@val) && (altv@val > retv@val))) {
										retv <- altv
									} 
								} else {
									if (is.na(retv@val) || is.nan(retv@val)) {
										# noop
									} else if (is.na(altv@val) || is.nan(altv@val) || (altv@val > retv@val)) {
										retv <- altv
									} 
								}
								vtag(retv) <- vtag
							} else {
								vtag <- paste0('max(',retv@vtag,', numeric, na.rm=',na.rm,')')

								if (na.rm) {
									if (is.na(retv@val) || is.nan(retv@val) || (!is.na(altv) && !is.nan(altv) && (altv > retv@val))) {
										val <- altv
										dvdx <- matrix(0,nrow=1,ncol=ncol(retv@dvdx))
									} 
								} else {
									if (is.na(retv@val) || is.nan(retv@val)) {
										# noop
									} else if (is.na(altv) || is.nan(altv) || (altv > retv@val)) {
										val <- altv
										dvdx <- matrix(0,nrow=1,ncol=ncol(retv@dvdx))
									} 
								}

								retv <- new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
								}
							}
						
						retv
					})

#' @rdname max
#' @aliases min
#' @aliases min,madness-class
setMethod("min", signature(x="madness"),
					function(x, ..., na.rm = FALSE) {
						xtag <- x@xtag
						val <- min(x@val,na.rm=na.rm)

						if (na.rm) { 
							dvdx <- x@dvdx[which.min(x@val),,drop=FALSE]
						} else {
							isok <- !(is.na(x@val) | is.nan(x@val))
							if (all(isok)) {
								dvdx <- x@dvdx[which.min(x@val),,drop=FALSE]
							} else {
								dvdx <- matrix(NA,nrow=1,ncol=ncol(x@dvdx))
							}
						}
						dvdx <- matrix(dvdx,nrow=1)

						vtag <- paste0('min(',x@vtag,', na.rm=',na.rm,')')
						varx <- x@varx

						retv <- new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
						if (length(list(...)) > 0) {
							altv <- min(...,na.rm=na.rm)

							# combine them:
							if (is(altv,'madness')) {
								xtag <- .check_common_xtag(retv,altv)
								varx <- .get_a_varx(retv,altv)
								vtag <- paste0('min(',retv@vtag,', ',altv@vtag,', na.rm=',na.rm,')')
								if (na.rm) {
									if (is.na(retv@val) || is.nan(retv@val) || (!is.na(altv@val) && !is.nan(altv@val) && (altv@val < retv@val))) {
										retv <- altv
									} 
								} else {
									if (is.na(retv@val) || is.nan(retv@val)) {
										# noop
									} else if (is.na(altv@val) || is.nan(altv@val) || (altv@val < retv@val)) {
										retv <- altv
									} 
								}
								vtag(retv) <- vtag
							} else {
								vtag <- paste0('min(',retv@vtag,', numeric, na.rm=',na.rm,')')

								if (na.rm) {
									if (is.na(retv@val) || is.nan(retv@val) || (!is.na(altv) && !is.nan(altv) && (altv < retv@val))) {
										val <- altv
										dvdx <- matrix(0,nrow=1,ncol=ncol(retv@dvdx))
									} 
								} else {
									if (is.na(retv@val) || is.nan(retv@val)) {
										# noop
									} else if (is.na(altv) || is.nan(altv) || (altv < retv@val)) {
										val <- altv
										dvdx <- matrix(0,nrow=1,ncol=ncol(retv@dvdx))
									} 
								}

								retv <- new("madness", val=val, dvdx=dvdx, vtag=vtag, xtag=xtag, varx=varx)
								}
							}
						
						retv
					})

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
