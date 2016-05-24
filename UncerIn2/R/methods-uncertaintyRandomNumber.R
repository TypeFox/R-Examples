## This file is part of the UncertaintyInterpolation 2.0 package.
##
## Copyright 2015 Tomas Burian


#' @title
#' Creates S4 object class UncertainPoints
#'
#' @description
#' Builds an uncertainty model based on the imprecision made by 
#' random number (from defined interval) over the input data. 
#' Input data must be type of S4 class \code{Points}. Output 
#' object is type of S4 class \code{UncertainPoints}.
#'
#' @param data Input data. S4 class of \code{Points}.
#' @param min The minimum value of the interval.
#' @param max The maximum value of the interval.
#' 
#' @return Returns an object of class \code{UncertainPoints}.
#' 
#' @seealso \code{\link[UncerIn2]{Points-class}}, \code{\link[UncerIn2]{UncertainPoints-class}}, \code{\link[UncerIn2]{uncertaintyInterpolation2-package}}
#'
#' @export
#' @importFrom stats runif


uncertaintyRandomNumber <- function (data, min=0, max=1)
{
  
  if(!(inherits(data, "Points"))){
    stop( paste("Dataset: ", deparse(substitute(testData))," is not S4 class of type Points." , sep=""))
  }
  
  randomNumber = runif(length(data@z), min, max)
  
  uncertaintyLower = data@z - randomNumber
  uncertaintyUpper = data@z + randomNumber 
  
  new("UncertainPoints", x=data@x, y=data@y, uncertaintyLower=uncertaintyLower, modalValue=data@z, uncertaintyUpper=uncertaintyUpper)
}