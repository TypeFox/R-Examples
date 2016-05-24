# /usr/bin/r
#
# Created: 2015.11.23
# Copyright: Steven E. Pav, 2015
# Author: Steven E. Pav <shabbychef@gmail.com>
# Comments: Steven E. Pav

#' @include AllClass.r
#' @include utils.r
NULL

#' @title Calculate Variance-Covariance Matrix for a model.
#'
#' @description 
#'
#' Returns the variance-covariance matrix of the parameters
#' computed by a \code{madness} object.
#'
#' @details
#'
#' Let \eqn{X}{X} represent some quantity which is estimated from
#' data. Let \eqn{\Sigma}{Sigma} be the (known or estimated)
#' variance-covariance matrix of \eqn{X}{X}. If \eqn{Y}{Y}
#' is some computed function of \eqn{X}{X}, then, by the 
#' Delta method (which is a first order Taylor approximation),
#' the variance-covariance matrix of \eqn{Y}{Y} is approximately
#' \deqn{\frac{\mathrm{d}Y}{\mathrm{d}{X}} \Sigma \left(\frac{\mathrm{d}Y}{\mathrm{d}{X}}\right)^{\top},}{(dY/dX) Sigma (dY/dX)',}
#' where the derivatives are defined over the 'unrolled' (or vectorized)
#' \eqn{Y}{Y} and \eqn{X}{X}. 
#'
#' Note that \eqn{Y}{Y} can represent a multidimensional quantity. Its
#' variance covariance matrix, however, is two dimensional, as it too
#' is defined over the 'unrolled' \eqn{Y}{Y}.
#'
#' @param object a \code{madness} object. A \code{varx} matrix must have
#' been set on the object, otherwise an error will be thrown.
#' @param ... additional arguments for method functions. Ignored here.
#' @export 
#' @method vcov madness
#' @return A matrix of the estimated covariances between the values being
#' estimated by the \code{madness} object. While \eqn{Y}{Y} may be 
#' multidimensional, the return value is a square matrix whose side length
#' is the number of elements of \eqn{Y}{Y}
#' @seealso \code{\link[stats]{vcov}}.
#' @examples 
#' y <- array(rnorm(2*3),dim=c(2,3))
#' dy <- matrix(rnorm(length(y)*2),ncol=2)
#' dx <- crossprod(matrix(rnorm(ncol(dy)*100),nrow=100))
#' obj <- madness(val=y,dvdx=dy,varx=dx)
#' print(vcov(obj))
#'
#' @template etc
vcov.madness <- function(object,...) {
	stopifnot(dim(object@varx)[1] != 0)
	retv <- object@dvdx %*% (object@varx %*% t(object@dvdx))
	invisible(retv)
}

#for vim modeline: (do not edit)
# vim:fdm=marker:fmr=FOLDUP,UNFOLD:cms=#%s:syn=r:ft=r
