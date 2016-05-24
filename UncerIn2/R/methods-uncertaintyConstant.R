## This file is part of the UncertaintyInterpolation 2.0 package.
##
## Copyright 2015 Tomas Burian


#' @title
#' Creates S4 object class UncertainPoints
#'
#' @description
#' Builds an uncertainty model based on the constant error over the input data. 
#' Input data must be type of S4 class \code{Points}. Output object 
#' is type of S4 class \code{UncertainPoints}.
#'
#' @param data Input data. S4 class of \code{Points}.
#' @param value Constant value of imprecision that defines the uncertainty.
#' 
#' @return Returns an object of class \code{UncertainPoints}.
#'
#' @seealso \code{\link[UncerIn2]{Points-class}}, \code{\link[UncerIn2]{UncertainPoints-class}}, \code{\link[UncerIn2]{uncertaintyInterpolation2-package}}
#' @export


uncertaintyConstant <- function(data, value = 2)
{
  
  if(!(inherits(data, "Points"))){
    stop( paste("Dataset: ", deparse(substitute(testData))," is not S4 class of type Points." , sep=""))
  }
  
  uncertaintyLower = data@z - value
  uncertaintyUpper = data@z + value  
  
  new("UncertainPoints", x=data@x, y=data@y, uncertaintyLower=uncertaintyLower, modalValue=data@z, uncertaintyUpper=uncertaintyUpper)
}
