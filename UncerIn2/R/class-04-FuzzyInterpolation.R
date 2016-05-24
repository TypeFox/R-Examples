## This file is part of the UncertaintyInterpolation 2.0 package.
##
## Copyright 2015 Tomas Burian


#' @title
#' S4 class Representing a FuzzyInterpolation
#' 
#' @description
#' Definition of S4 class for output data from fuzzy (kriging) 
#' interpolation function of package \code{UncerIn2}. Where 
#' \code{x, y} should represent coordinates and \code{minimal, 
#' modalValue, maximal} the values of uncertainty.
#'
#' @section Slots:
#'  \describe{
#'    \item{\code{x}:}{Input numeric data (number specyfing the x grid coordinate).}
#'    \item{\code{y}:}{Input numeric data (number specyfing the y grid coordinate).}
#'    \item{\code{minimal}:}{Defined model values of uncertaintyLower - the bottom calculated part.}
#'    \item{\code{modalValue}:}{Defined model values of modalValue - the mean values.}
#'    \item{\code{maximal}:}{Defined model values of uncertaintyUpper - the upper calculated part.}
#'  }
#'
#' @exportClass FuzzyInterpolation
#' @name FuzzyInterpolation-class
#' @rdname FuzzyInterpolation-class
#' @seealso \code{\link[UncerIn2]{uncertaintyInterpolation2-package}}
#' @docType class
#' @family FuzzyInterpolation-method
#' @examples
#' showClass("FuzzyInterpolation")


setClass(Class = "FuzzyInterpolation", 
         representation=representation(
           x = "numeric",
           y = "numeric",
           minimal = "numeric",
           modalValue = "numeric",
           maximal = "numeric"
         ),
         
         prototype=prototype(
           x = c(),
           y = c(),
           minimal = c(),
           modalValue = c(),
           maximal = c()
         ),
         
         validity=function(object){
           if(length(object@x) != length(object@y) || length(object@x) != length(object@modalValue) || 
                length(object@x) != length(object@minimal) || length(object@x) != length(object@maximal)){
             return("Length of all parameters (x,y,z) should be equal.")
           }
           
           # OK
           return(TRUE)
         }         
)