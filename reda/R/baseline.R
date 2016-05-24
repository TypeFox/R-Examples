################################################################################
##
##   R package reda by Wenjie Wang, Haoda Fu, and Jun Yan
##   Copyright (C) 2015
##
##   This file is part of the R package reda.
##
##   The R package reda is free software: You can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   any later version (at your option). See the GNU General Public License
##   at <http://www.gnu.org/licenses/> for details.
##
##   The R package reda is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##
################################################################################


## collation after class.R
#' @include class.R 
NULL


#' Estimated Coefficients of Baseline Rate Function
#' 
#' An S4 class generic function that returns the estimated coefficients
#' of baseline rate function. For \code{\link{rateReg-class}} object,
#' it returns either coefficients of pieceswise (including one piece)
#' constant rate function or coefficients of B-spline bases.
#' 
#' @param object An object used to dispatch a method.
#' @param ... Other arguments for future usage.
#' @return A named numeric vector.
#' @aliases baseRate,rateReg-method
#' @examples
#' ## See examples given in function rateReg.
#' @seealso
#' \code{\link{rateReg}} for model fitting;
#' \code{\link{summary,rateReg-method}} for summary of a fitted model.
#' @export
setGeneric(name = "baseRate",
           def = function(object, ...) {
               standardGeneric("baseRate")
           })


#' @describeIn baseRate Extract estiamted coefficients of
#' baseline rate function from \code{rateReg-class} object.
#' @export
setMethod(f = "baseRate", signature = "rateReg",
          definition = function(object, ...) {
              object@estimates$alpha[, 1]
          })

