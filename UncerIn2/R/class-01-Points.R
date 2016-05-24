## This file is part of the UncertaintyInterpolation 2.0 package.
##
## Copyright 2015 Tomas Burian


#' @title
#' S4 class Representing a class Points
#' 
#' @description
#' Definition of class for input data of numeric format \code{x/y/z}. 
#' Also defines input data format into the functions of package \code{UncerIn2}, 
#' which are building uncertainty models over input data. Where \code{x, y} should 
#' represent coordinates and \code{z} variable.
#'
#' @section Slots:
#'  \describe{
#'    \item{\code{x}:}{Input numeric data (number specyfing the x coordinate).}
#'    \item{\code{y}:}{Input numeric data (number specyfing the y coordinate).}
#'    \item{\code{z}:}{Input numeric data (number specyfing the values of variable).}
#'  }
#'
#' @exportClass Points
#' @import methods
#' @import Rcpp
#' @name Points-class
#' @rdname Points-class
#' @seealso \code{\link[UncerIn2]{uncertaintyInterpolation2-package}}
#' @docType class
#' @family Points-method
#' @examples
#' showClass("Points")


setClass(Class = "Points", 
         representation=representation(
           x = "numeric",
           y = "numeric",
           z = "numeric"
         ),
         
         prototype=prototype(
           x = c(),
           y = c(),
           z = c()
         ),
         
         validity=function(object){
           if(length(object@x) != length(object@y) || length(object@x) != length(object@z)){
             return("Length of all parameters (x,y,z) should be equal.")
           }
           
           # OK
           return(TRUE)
         }
)


#' @title
#' Creates S4 object class Points
#'
#' @description
#' This function creates an object of S4 class \code{Points}.
#'
#' @param x Input numeric data (number specyfing the x coordinate).
#' @param y Input numeric data (number specyfing the y coordinate).
#' @param z Input numeric data (number specyfing the values of variable).
#' 
#' @return Returns an object of class \code{Points}.
#'
#' @export


Points <- function(x, y, z)
{
  new("Points", x = x, y = y, z = z)
}