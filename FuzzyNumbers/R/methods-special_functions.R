## This file is part of the FuzzyNumbers library.
##
## Copyright 2012-2014 Marek Gagolewski
##
##
## FuzzyNumbers is free software: you can redistribute it and/or modify
## it under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## FuzzyNumbers is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU Lesser General Public License for more details.
##
## You should have received a copy of the GNU Lesser General Public License
## along with FuzzyNumbers. If not, see <http://www.gnu.org/licenses/>.



#' @title
#' Arc-tangent
#'
#' @description
#' The arc-tangent of two arguments arctan2(y, x) returns the angle between 
#' the x-axis and the vector from the origin to (x, y) for PiecewiseLinearFuzzyNumbers.
#' 
#' @details
#' Note that resulting values are no longer from interval [-pi,pi] but 
#' [-1.5pi,pi], in order to provide valid fuzzy numbers as result.
#' 
#' 
#' @param y a PiecewiseLinearFuzzyNumber
#' @param x a PiecewiseLinearFuzzyNumber
#'
#' @return Returns a fuzzy number of the class \linkS4class{PiecewiseLinearFuzzyNumber}
#' indicating the angle specified by the input fuzzy numbers. The range of results is
#' [-1.5pi,pi].  
#'
#' @exportMethod arctan2
#' @docType methods
#' @name arctan2
#' @family special-method
#' @family PiecewiseLinearFuzzyNumber-method
#' @aliases arctan2,PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber-method
#' 
#' @usage
#' \S4method{arctan2}{PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber}(y, x)
#' 
#' @examples 
#' y = as.PiecewiseLinearFuzzyNumber(TriangularFuzzyNumber(-2, 3, 5), knot.n = 9)
#' x = as.PiecewiseLinearFuzzyNumber(TriangularFuzzyNumber(-4.8, -4, 1.5), knot.n = 9)
#' arctan2(y,x)
setGeneric("arctan2",
           function(y, x) standardGeneric("arctan2"))
setMethod(
  f="arctan2",
  signature(y="PiecewiseLinearFuzzyNumber", x="PiecewiseLinearFuzzyNumber"),
  definition=function(y, x)
  {
    fns = fuzzyNumberUnification(y,x)
    y = fns[[1]]
    x = fns[[2]]
    
    knot.alpha = (y@knot.alpha)

    yl <- c(y@a1,y@knot.left,y@a2)
    xl <- c(x@a1,x@knot.left,x@a2)
    yr <- rev(c(y@a3,y@knot.right,y@a4))
    xr <- rev(c(x@a3,x@knot.right,x@a4))
    
    modified = FALSE
    mod <- ((yl<=0 & 0<=yr) & (xr<=0))
    #if there are no zeroes than there is no need for modified atan2
    modified = (length((which(mod == TRUE)))!=0)

    #identify alpha cuts that both contain zero
    zero <- ((yl<=0 & 0<=yr) & (xl<=0 & 0<=xr))

    if(modified){
      p1 <- apply(cbind(yl,xl), 1, function(x){arctan2m(x[1],x[2])})
      p2 <- apply(cbind(yl,xr), 1, function(x){arctan2m(x[1],x[2])})
      p3 <- apply(cbind(yr,xl), 1, function(x){arctan2m(x[1],x[2])})
      p4 <- apply(cbind(yr,xr), 1, function(x){arctan2m(x[1],x[2])})
      
      #for intervals with zeroes set fixed limiting values
      p1[zero] <- -1.5*pi 
      p4[zero] <-  0.5*pi
    }else{
      p1 <- apply(cbind(yl,xl), 1, function(x){arctan2e(x[1],x[2])})
      p2 <- apply(cbind(yl,xr), 1, function(x){arctan2e(x[1],x[2])})
      p3 <- apply(cbind(yr,xl), 1, function(x){arctan2e(x[1],x[2])})
      p4 <- apply(cbind(yr,xr), 1, function(x){arctan2e(x[1],x[2])})
      
      #for intervals with zeroes set fixed limiting values
      p1[zero] <- -pi 
      p4[zero] <-  pi
    }
    
    result <-  PiecewiseLinearFuzzyNumber(knot.alpha=knot.alpha,
                                          knot.left=pmin(p1, p2, p3, p4),
                                          knot.right=rev(pmax(p1, p2, p3, p4)))
    
    return(result)
  }
)


#' @title
#' Integer power of fuzzy number
#'
#' @description
#' For fuzzy numbers the equality of \code{X*X == X^2} does not hold. 
#' 
#' @details
#' This function calculates integer power of a PiecewiseLinearFuzzyNumber according
#' to the reference below.
#' 
#' 
#' @param e1 a PiecewiseLinearFuzzyNumber
#' @param e2 numeric (if it is not integer it will be converted by function
#' as.integer())
#'
#' @return Returns a fuzzy number of the class \linkS4class{PiecewiseLinearFuzzyNumber}
#' indicating \code{e1^e2}. 
#'
#' @references KAUFMANN, A., GUPTA, M. M. (1985) Introduction to Fuzzy Arithmetic. 
#' New York : Van Nostrand Reinhold Company. ISBN 044230079. 
#'
#' @docType methods
#' @family extension_principle
#' @family PiecewiseLinearFuzzyNumber-method
#' @aliases ^,PiecewiseLinearFuzzyNumber,numeric-method
#' 
#' @usage
#' \S4method{^}{PiecewiseLinearFuzzyNumber,numeric}(e1, e2)
#' 
#' @examples 
#' x = as.PiecewiseLinearFuzzyNumber(TriangularFuzzyNumber(-2, 1, 9), knot.n = 2)
#' x^2
#' x^3
setMethod(
  "^",
  signature(e1="PiecewiseLinearFuzzyNumber", e2 ="numeric"),
  function(e1, e2)
  {
    e2 = as.integer(e2)
    knot.alpha = (e1@knot.alpha)
    
    xl <- (c(e1@a1,e1@knot.left,e1@a2))
    xr <- (rev(c(e1@a3,e1@knot.right,e1@a4)))
    
    p1 <- xl^e2
    p2 <- xr^e2
    
    if(e2%%2 == 0){
      
      zero <- (xl<=0 & 0<=xr)
      
      min <- pmin(p1, p2)
      min[zero] <- 0
      max <- rev(pmax(p1, p2))
      
    }else{
      min <- pmin(p1, p2)
      max <- rev(pmax(p1, p2))
    }
    
    result <-  PiecewiseLinearFuzzyNumber(knot.alpha=knot.alpha,
                                          knot.left=min,
                                          knot.right=max)
    
    return(result)
  }
)




#following methods are merely helpers to correctly implement arctan2 for PiecewiseLinearFuzzyNumbers
setGeneric("arctan2e",
           function(y, x) standardGeneric("arctan2e"))
setMethod(
  f="arctan2e",
  signature(y="numeric", x="numeric"),
  definition=function(y, x)
  {
    result = 0;
    
    if(x>0){
      result = atan(y/x);
    }
    else if(x<0 & 0<=y){
      result = atan(y/x) + pi;
    }
    else if(x<0 & 0>y){
      result = atan(y/x) - pi;
    }
    else if(x==0 & y>0){
      result = pi/2;
    }
    else if(x==0 & y<0){
      result = - pi/2;
    }
    else if(x==0 & y==0){
      result = 0;
    }
    
    return (result)
  }
)

setGeneric("arctan2m",
           function(y, x) standardGeneric("arctan2m"))
setMethod(
  f="arctan2m",
  signature(y="numeric", x="numeric"),
  definition=function(y, x)
  {
    result = 0;
    
    if(x>0){
      result = atan(y/x);
    }
    else if(x<0 & 0<=y){
      result = atan(y/x) - pi;
    }
    else if(x<0 & 0>y){
      result = atan(y/x) - pi;
    }
    else if(x==0 & y>0){
      result = pi/2;
    }
    else if(x==0 & y<0){
      result = - pi/2;
    }
    else if(x==0 & y==0){
      result = 0;
    }
    
    return (result)
  }
)
