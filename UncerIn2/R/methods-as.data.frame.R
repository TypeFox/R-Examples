## This file is part of the UncertaintyInterpolation 2.0 package.
##
## Copyright 2015 Tomas Burian


#' @title
#' Transformation of S4 class Points or UncertainPoints into data.frame
#'
#' @description
#' This function provides the transformation of S4 object class \code{UncertainPoints} or 
#' \code{Points} into the \code{data.frame} data format. 
#'
#' @param data Input data type of S4 object of class UncertainPoints or Points.
#' 
#' @usage
#' \S4method{as.dataframe}{UncertainPoints}(data)
#'
#' \S4method{as.dataframe}{Points}(data)
#'
#' @return Returns an object of class \code{data.frame}.
#'
#' @name as.dataframe
#' @docType methods
#' @rdname as.dataframe
#' @seealso \code{\link[UncerIn2]{Points-class}}, 
#' \code{\link[UncerIn2]{Points}}, \code{\link[UncerIn2]{UncertainPoints-class}}, 
#' \code{\link[UncerIn2]{uncertaintyInterpolation2-package}}
#' @aliases as.dataframe,UncertainPoints-method
#'          as.dataframe,Points-method
#' @exportMethod as.dataframe


setGeneric("as.dataframe",
            function(data) 
            standardGeneric("as.dataframe")
           )
	

setMethod(f="as.dataframe",
          signature(data="UncertainPoints"),
          definition=function(data)
        {
        
        a <- cbind(data@x, data@y, data@uncertaintyLower, data@modalValue, data@uncertaintyUpper) 
        b = as.data.frame(a)
        
        colnames(b) <-  c("x", "y", "uncertaintyLower", "modalValue", "uncertaintyUpper")           
        result = b
        
        return(result)
        }
)


setMethod(f="as.dataframe",
          signature(data="Points"),
          definition=function(data)
          {
            
            a <- cbind(data@x, data@y, data@z) 
            b = as.data.frame(a)
            
            colnames(b) <-  c("x", "y", "z")           
            result = b
            
            return(result)
          }
)