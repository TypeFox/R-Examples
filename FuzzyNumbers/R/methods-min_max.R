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
#' Maximum of fuzzy nubmers
#'
#' @description
#' Determines maximum fuzzy number based on two inputs. 
#' 
#' @details
#' The resulting 
#' PiecewiseLinearFuzzyNumber is only an approximation of the result
#' it might not be very precise for small number of knots (see examples).
#' 
#' @references KAUFMANN, A., GUPTA, M. M. (1985) Introduction to Fuzzy Arithmetic.
#'  New York : Van Nostrand Reinhold Company. ISBN 044230079. 
#' 
#' @param e1 a PiecewiseLinearFuzzyNumber
#' @param e2 a PiecewiseLinearFuzzyNumber
#'
#' @return Returns a PiecewiseLinearFuzzyNumber representing maximal
#' value of the inputs.
#'
#' @exportMethod maximum
#' @docType methods
#' @name maximum
#' @family min_max-operators
#' @family PiecewiseLinearFuzzyNumber-method
#' @aliases maximum,PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber-method
#' 
#' @usage
#' \S4method{maximum}{PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber}(e1, e2)
#' 
#' @examples
#' # example with low number of knots, showing the approximate nature
#' # of the result
#' x = as.PiecewiseLinearFuzzyNumber(TriangularFuzzyNumber(-4.8, -3 , -1.5))
#' y = as.PiecewiseLinearFuzzyNumber(TriangularFuzzyNumber(-5.5, -2.5, -1.1))
#' maxFN = maximum(x,y)
#' min = min(alphacut(x,0)[1],alphacut(y,0)[1])
#' max = max(alphacut(x,0)[2],alphacut(y,0)[2])
#' plot(x, col="red", xlim=c(min,max))
#' plot(y, col="blue", add=TRUE)
#' plot(maxFN, col="green", add=TRUE)
#' 
#' # example with high number of knots, that does not suffer 
#' # from the approximate nature of the result
#' x = as.PiecewiseLinearFuzzyNumber(TriangularFuzzyNumber(-4.8, -3 , -1.5), knot.n = 9)
#' y = as.PiecewiseLinearFuzzyNumber(TriangularFuzzyNumber(-5.5, -2.5, -1.1), knot.n = 9)
#' maxFN = maximum(x,y)
#' min = min(alphacut(x,0)[1],alphacut(y,0)[1])
#' max = max(alphacut(x,0)[2],alphacut(y,0)[2])
#' plot(x, col="red", xlim=c(min,max))
#' plot(y, col="blue", add=TRUE)
#' plot(maxFN, col="green", add=TRUE)
setGeneric("maximum",
           function(e1,e2) standardGeneric("maximum"))

setMethod(
  "maximum",
  signature(e1 = "PiecewiseLinearFuzzyNumber", e2 = "PiecewiseLinearFuzzyNumber"),
  function (e1, e2)
  {
    fns = fuzzyNumberUnification(e1,e2)
    e1 = fns[[1]]
    e2 = fns[[2]]
    
    e1l <- c(e1@a1,e1@knot.left,e1@a2)
    e2l <- c(e2@a1,e2@knot.left,e2@a2)
    e1r <- rev(c(e1@a3,e1@knot.right,e1@a4))
    e2r <- rev(c(e2@a3,e2@knot.right,e2@a4))
    
    knot.alpha <- unique(sort(c(e1@knot.alpha, e2@knot.alpha)))
    
    result <- PiecewiseLinearFuzzyNumber(knot.alpha=knot.alpha,
                               knot.left=pmax(e1l,e2l),
                               knot.right=rev(pmax(e1r,e2r)))
    return(result)
  }
)


#' @title
#' Minimum of fuzzy nubmers
#'
#' @description
#' Determines minimum fuzzy number based on two inputs.
#' 
#' @details
#' The resulting 
#' PiecewiseLinearFuzzyNumber is only an approximation of the result
#' it might not be very precise for small number of knots (see examples).
#' 
#' @references KAUFMANN, A., GUPTA, M. M. (1985) Introduction to Fuzzy Arithmetic.
#' New York : Van Nostrand Reinhold Company. ISBN 044230079. 
#' 
#' @param e1 a PiecewiseLinearFuzzyNumber
#' @param e2 a PiecewiseLinearFuzzyNumber
#'
#' @return Returns a PiecewiseLinearFuzzyNumber representing maximal
#' value of the inputs.
#'
#' @exportMethod minimum
#' @docType methods
#' @name minimum
#' @family min_max-operators
#' @family PiecewiseLinearFuzzyNumber-method
#' @aliases minimum,PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber-method
#' 
#' @usage
#' \S4method{minimum}{PiecewiseLinearFuzzyNumber,PiecewiseLinearFuzzyNumber}(e1, e2)
#' 
#' @examples
#' # example with low number of knots, showing the approximate nature
#' # of the result
#' x = as.PiecewiseLinearFuzzyNumber(TriangularFuzzyNumber(-4.8, -3 , -1.5))
#' y = as.PiecewiseLinearFuzzyNumber(TriangularFuzzyNumber(-5.5, -2.5, -1.1))
#' minFN = minimum(x,y)
#' min = min(alphacut(x,0)[1],alphacut(y,0)[1])
#' max = max(alphacut(x,0)[2],alphacut(y,0)[2])
#' plot(x, col="red", xlim=c(min,max))
#' plot(y, col="blue", add=TRUE)
#' plot(minFN, col="green", add=TRUE)
#' 
#' # example with high number of knots, that does not suffer 
#' # from the approximate nature of the result
#' x = as.PiecewiseLinearFuzzyNumber(TriangularFuzzyNumber(-4.8, -3 , -1.5), knot.n = 9)
#' y = as.PiecewiseLinearFuzzyNumber(TriangularFuzzyNumber(-5.5, -2.5, -1.1), knot.n = 9)
#' minFN = minimum(x,y)
#' min = min(alphacut(x,0)[1],alphacut(y,0)[1])
#' max = max(alphacut(x,0)[2],alphacut(y,0)[2])
#' plot(x, col="red", xlim=c(min,max))
#' plot(y, col="blue", add=TRUE)
#' plot(minFN, col="green", add=TRUE)
setGeneric("minimum",
           function(e1,e2) standardGeneric("minimum"))

setMethod(
  "minimum",
  signature(e1 = "PiecewiseLinearFuzzyNumber", e2 = "PiecewiseLinearFuzzyNumber"),
  function (e1, e2)
  {
    fns = fuzzyNumberUnification(e1,e2)
    e1 = fns[[1]]
    e2 = fns[[2]]
    
    e1l <- c(e1@a1,e1@knot.left,e1@a2)
    e2l <- c(e2@a1,e2@knot.left,e2@a2)
    e1r <- rev(c(e1@a3,e1@knot.right,e1@a4))
    e2r <- rev(c(e2@a3,e2@knot.right,e2@a4))
    
    knot.alpha <- unique(sort(c(e1@knot.alpha, e2@knot.alpha)))
    
    result <- PiecewiseLinearFuzzyNumber(knot.alpha=knot.alpha,
                                         knot.left=pmin(e1l,e2l),
                                         knot.right=rev(pmin(e1r,e2r)))
    return(result)
  }
)
