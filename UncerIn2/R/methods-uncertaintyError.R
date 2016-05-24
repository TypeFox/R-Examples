## This file is part of the UncertaintyInterpolation 2.0 package.
##
## Copyright 2015 Tomas Burian


#' @title
#' Creates S4 object class UncertainPoints
#'
#' @description
#' Builds an uncertainty model over the input data based on the spatially 
#' correlated errors. Input data must be type of S4 class \code{Points}. 
#' Output object is type of S4 class \code{UncertainPoints}.
#' 
#' @param input Input data. S4 class of \code{Points}
#' @param T optional vector of time coordinates, T must always be an equidistant vector. Instead of T=seq(from=From, by=By, len=Len) one may also write T=c(From, By, Len).
#' @param grid Logical; RandomFields can find itself the correct value in nearly all cases, so that usually grid need not be given.
#' @param distances Another alternative to pass the (relative) coordinates.
#' @param dim Only used if distances are given.
#' @param data For conditional simulation and random imputing only. If data is missing, unconditional simulation is performed.
#' @param given Optional, matrix or list. If given matrix then the coordinates can be given separately, namely by given where, in each row, a single location is given.
#' 
#' @return Returns an object of class \code{UncertainPoints}.
#'
#' @details For the calculations of spatially correlated errors was used package RandomFields.
#'
#' @seealso \code{\link[UncerIn2]{Points-class}}, \code{\link[UncerIn2]{UncertainPoints-class}}, \code{\link[RandomFields]{RFsimulate}}, \code{\link[UncerIn2]{uncertaintyInterpolation2-package}}
#' @import RandomFields
#' @export


uncertaintyError <- function(input, T=NULL, grid=FALSE, distances=NULL, dim=NULL, data=NULL, given=NULL)   
{
  
  if(!(inherits(input, "Points"))){
    stop( paste("Dataset: ", deparse(substitute(testData))," is not S4 class of type Points." , sep=""))
  }
  
  model <- RMstable(alpha=1.9, scale=0.4)
  error <- suppressMessages(RFsimulate(model=model, x=input@x, y=input@y, T=T, grid=grid, 
                                       distances=distances, dim=dim, data=data, given=given, n=1))
  
  var = abs(error@data[,1])
  
  uncertaintyLower = input@z - var
  uncertaintyUpper = input@z + var
  
  new("UncertainPoints", x=input@x, y=input@y, uncertaintyLower=uncertaintyLower, modalValue=input@z, uncertaintyUpper=uncertaintyUpper)
}