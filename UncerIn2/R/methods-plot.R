## This file is part of the UncertaintyInterpolation 2.0 package.
##
## Copyright 2015 Tomas Burian


#' @title
#' Plotting S4 class UncertainInterpolation
#'
#' @description
#' This function provides the plotting of S4 object class \code{UncertainInterpolation}. 
#'
#' @param object Input data type of S4 object class UncertainInterpolation.
#' @param attr1 First plotting atrribute.
#' @param attr2 Second plotting atrribute.
#' @param attr3 Third plotting atrribute.
#' @param cuts Number of cuts.
#' @param pretty Logical value \code{TRUE/FALSE.}(choose colour breaks at pretty numbers?)
#' 
#' @usage
#' \S4method{Plot}{UncertainInterpolation}(object, attr1, attr2, attr3, cuts, pretty)
#' 
#' \S4method{Plot}{UncertainInterpolation}(object, attr1 = "uncertaintyLower", attr2 = "modalValue", 
#' attr3 = "uncertaintyUpper", cuts = 10, pretty=TRUE)
#' 
#' @seealso \code{\link[UncerIn2]{UncertainInterpolation-class}}, \code{\link[UncerIn2]{uncertaintyInterpolation2-package}}
#' 
#' @name Plot
#' @docType methods
#' @rdname Plot
#' @aliases Plot,UncertainInterpolation-method
#' 
#' @exportMethod Plot


setGeneric("Plot",
           function(object, ...)
           standardGeneric("Plot")
)


setMethod("Plot",
          signature(object = "UncertainInterpolation"),
          definition = function(object, attr1 = "uncertaintyLower", attr2 = "modalValue", attr3 = "uncertaintyUpper"
                                , cuts = 10, pretty=TRUE)           
          {
            
            a = as.UncertainPoints(object)
            a = as.dataframe(a)
            
            gridded(a)=~x+y         
          
            spplot(a, c(attr1, attr2, attr3), names.attr= c(attr1, attr2, attr3),
                   colorkey=list(space="bottom"), layout=c(3,1), cuts = cuts, pretty=pretty)            
          }
)