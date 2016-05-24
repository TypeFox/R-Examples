## This file is part of the UncertaintyInterpolation 2.0 package.
##
## Copyright 2015 Tomas Burian


#' @title
#' Transformation of S4 class UncertainInterpolation into UncertainPoints
#'
#' @description
#' This function provides the transformation of S4 object class \code{UncertainInterpolation} 
#' into the S4 class \code{UncertainPoints}.. 
#' 
#' @param object Input data type of S4 object class UncertainInterpolation.
#' @param grid Input grid (default set FALSE)
#'
#' @usage
#' \S4method{as.UncertainPoints}{UncertainInterpolation}(object,grid)
#'
#' @return Returns an object of class \code{UncertainPoints}.
#' 
#' @seealso \code{\link[UncerIn2]{UncertainInterpolation-class}}, \code{\link[UncerIn2]{UncertainPoints-class}}, \code{\link[UncerIn2]{uncertaintyInterpolation2-package}}
#' 
#' @name as.UncertainPoints
#' @docType methods
#' @rdname as.UncertainPoints
#' @aliases as.UncertainPoints,UncertainInterpolation-method
#' @exportMethod as.UncertainPoints


setGeneric("as.UncertainPoints",
            function(object, ...) 
            standardGeneric("as.UncertainPoints")
           )

	
setMethod(f="as.UncertainPoints",
          signature(object="UncertainInterpolation"),
          definition=function(object, grid = FALSE)
        {
 
              a = new("UncertainPoints", x=object@x, y=object@y, uncertaintyLower=object@uncertaintyLower, modalValue=object@modalValue, uncertaintyUpper=object@uncertaintyUpper)
              
              return(a)    
        }
)