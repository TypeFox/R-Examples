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
# Created: 2016.01.14
# Copyright: Steven E. Pav, 2016
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

#' @include AllClass.r
#' @include utils.r
NULL

#' @title Convert a madness object into an objective value with gradient
#'
#' @description 
#'
#' Given a \code{madness} object representing a scalar value, strip out
#' that value and attach an attribute of its derivative as a gradient.
#' This is a convenience method that simplifies construction of objective
#' functions for optimization routines.
#'
#' @usage
#'
#' to_objective(X)
#'
#' @param X a \code{madness} object representing a scalar. 
#' @note An error will be thrown if the value is not a scalar.
#' @return A scalar numeric with a \code{gradient} attribute of the derivative.
#' @template etc
#' @examples 
#' # an objective function for matrix factorization with penalty:
#' fitfun <- function(R,L,Y,nu=-0.1) {
#'	Rmad <- madness(R)
#'	dim(Rmad) <- c(ncol(L),ncol(Y))
#'	Err <- Y - L %*% Rmad
#'	penalty <- sum(exp(nu * Rmad))
#'	fit <- norm(Err,'f') + penalty
#'	to_objective(fit)
#'	}
#' set.seed(1234)
#' L <- array(runif(30*5),dim=c(30,5)) 
#' Y <- array(runif(nrow(L)*20),dim=c(nrow(L),20))
#' R0 <- array(runif(ncol(L)*ncol(Y)),dim=c(ncol(L),ncol(Y)))
#' obj0 <- fitfun(R0,L,Y)
#' fooz <- nlm(fitfun, R0, L, Y, iterlim=3)
#' @export
#' @rdname to_objective
to_objective <- function(X) {
	retv <- val(X)
	stopifnot(length(retv)==1)
	retv <- as.numeric(retv)
	attr(retv, "gradient") <- t(dvdx(X))
	retv
}

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
