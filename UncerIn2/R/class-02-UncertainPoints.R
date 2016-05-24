## This file is part of the UncertaintyInterpolation 2.0 package.
##
## Copyright 2015 Tomas Burian


#' @title
#' S4 class Representing a UncertainPoints
#' 
#' @description
#' Definition of S4 class for output data from functions of package 
#' \code{UncerIn2}, which are building uncertainty models over input data. 
#' Also defines input data format into the interpolation function of package 
#' \code{UncerIn2}. Where \code{x, y} should represent coordinates and \code{uncertaintyLower, 
#' modalValue, uncertaintyUpper} the variables of uncertainty model.
#'
#' @section Slots:
#'  \describe{
#'    \item{\code{x}:}{Input numeric data (number specyfing the x coordinate).}
#'    \item{\code{y}:}{Input numeric data (number specyfing the y coordinate).}
#'    \item{\code{uncertaintyLower}:}{Defined model values of uncertaintyLower - the bottom calculated part.}
#'    \item{\code{modalValue}:}{Defined model values of modalValue - the mean values.}
#'    \item{\code{uncertaintyUpper}:}{Defined model values of uncertaintyUpper - the upper calculated part.}
#'  }
#'
#' @exportClass UncertainPoints
#' @name UncertainPoints-class
#' @rdname UncertainPoints-class
#' @seealso \code{\link[UncerIn2]{uncertaintyInterpolation2-package}}
#' @docType class
#' @family UncertainPoints-method
#' @examples
#' showClass("UncertainPoints")


setClass(Class = "UncertainPoints", 
         representation=representation(
           x = "numeric",
           y = "numeric",
           uncertaintyLower = "numeric",
           modalValue = "numeric",
           uncertaintyUpper = "numeric"
         ),
         
         prototype=prototype(
           x = c(),
           y = c(),
           uncertaintyLower = c(),
           modalValue = c(),
           uncertaintyUpper = c()
         ),
         
         validity=function(object){
           if(length(object@x) != length(object@y) || length(object@x) != length(object@modalValue) || 
              length(object@x) != length(object@uncertaintyLower) || length(object@x) != length(object@uncertaintyUpper)){
             return("Length of all parameters (x,y,z) should be equal.")
           }
           
           if(object@uncertaintyLower >= object@modalValue || object@modalValue >= object@uncertaintyUpper){
              return("Error in mathematical inequality of the values: should be uncertaintyLower <= modalValue <= uncertaintyUpper")
           }
           
          # OK
           return(TRUE)
         }         
)