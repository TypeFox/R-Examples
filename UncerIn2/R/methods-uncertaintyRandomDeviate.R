## This file is part of the UncertaintyInterpolation 2.0 package.
##
## Copyright 2015 Tomas Burian


#' @title
#' Creates S4 object class UncertainPoints
#'
#' @description
#' Builds an uncertainty model based on the imprecision made by random 
#' deviates over the input data. Input data must be type of S4 class 
#' \code{Points}. Output object is type of S4 class \code{UncertainPoints}.
#' 
#' @param data Input data. S4 class of \code{Points}.
#' @param mean The mean value of input data values.
#' @param sd The standart deviation of input data values.
#' 
#' @return Returns an object of class \code{UncertainPoints}.
#' 
#' @seealso \code{\link[UncerIn2]{Points-class}}, \code{\link[UncerIn2]{UncertainPoints-class}}, \code{\link[UncerIn2]{uncertaintyInterpolation2-package}}
#' 
#' @export
#' @importFrom stats rnorm


uncertaintyRandomDeviate <- function(data, mean=0, sd=1)
{
  
  if(!(inherits(data, "Points"))){
    stop( paste("Dataset: ", deparse(substitute(testData))," is not S4 class of type Points." , sep=""))
  }
  
  value = rnorm(length(data@z), mean, sd)
  modify = abs(value)

  uncertaintyLower = data@z - modify
  uncertaintyUpper = data@z + modify 
  
  new("UncertainPoints", x=data@x, y=data@y, uncertaintyLower=uncertaintyLower, modalValue=data@z, uncertaintyUpper=uncertaintyUpper)
}