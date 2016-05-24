## This file is part of the UncertaintyInterpolation 2.0 package.
##
## Copyright 2015 Tomas Burian


#' @title
#' Generation of grid with defined number of cels
#'
#' @description
#' This method creates a grid defined in two possible ways - the user should decide if spatial grid is 
#' required (\code{gridded = TRUE/FALSE}). The size of the grid should be also defined by the user.
#'
#' @param object Input data type of S4 object class \code{Points}.
#' @param gridded Logical value \code{TRUE/FALSE}. If true, then Spatial grid is made.
#' @param numberOfCellsX The number of cells on axle X.
#' @param numberOfCellsY The number of cells on axle Y.
#'
#' @usage
#' \S4method{Grid.def}{Points}(object,gridded,numberOfCellsX,numberOfCellsY)
#'
#' @return Returns an object of class \code{SpatialPixels} or \code{data.frame}, depends on whether 
#' the spatial grid is required.
#' 
#' @seealso \code{\link[UncerIn2]{Points-class}}, \code{\link[UncerIn2]{Grid.box}}, \code{\link[UncerIn2]{Grid.interpolation}}, \code{\link[UncerIn2]{uncertaintyInterpolation2-package}}
#' 
#' @name Grid.def
#' @docType methods
#' @rdname Grid.def
#' @aliases Grid.def,Points-method
#' @exportMethod Grid.def
#' @import geoR


setGeneric("Grid.def",
           function(object, ...)
           standardGeneric("Grid.def")
)


setMethod("Grid.def",
          signature(object = "Points"),
          definition = function(object, gridded = FALSE, numberOfCellsX = 100, numberOfCellsY = 100)
          {
            
            if (identical(gridded, FALSE)) 
              {
                minX = min(object@x)
                maxX = max(object@x)
                minY = min(object@y)
                maxY = max(object@y)
                          
                stepX = (maxX-minX)/numberOfCellsX
                stepY = (maxY-minY)/numberOfCellsY
                          
                grid = expand.grid(x = seq(minX, maxX, stepX), y = seq(minY, maxY, stepY))
                          
                result=grid
                return(result)
              }
            
            else if (identical(gridded, TRUE)) 
            {
              minX = min(object@x)
              maxX = max(object@x)
              minY = min(object@y)
              maxY = max(object@y)
              
              stepX = (maxX-minX)/numberOfCellsX
              stepY = (maxY-minY)/numberOfCellsY
              
              grid = expand.grid(x = seq(minX, maxX, stepX), y = seq(minY, maxY, stepY))
              
              names(grid) = c("x", "y")
              gridded(grid) = ~x+y
              
              result=grid
              return(result)   # class SpatialPixels 
            }
            
            else {
              stop("argument `gridded` should be a logical TRUE or FALSE")
              
            }
    }
)


#' @title
#' Generation of grid defined by surrounding box
#'
#' @description
#' This method creates a grid defined in two possible ways - the user should decide if spatial grid is 
#' required (\code{gridded = TRUE/FALSE}). The size of the grid is recognized by the input data coordinates. 
#' User should define size of the single cell.
#'
#' @param object Input data type of S4 object class \code{Points}.
#' @param gridded Logical value \code{TRUE/FALSE}. If true, then Spatial grid is made.
#' @param cellsize Numeric value for the size of cell.
#'
#' @usage
#' \S4method{Grid.box}{Points}(object,gridded,cellsize)
#' 
#' @return Returns an object of class \code{SpatialPixels} or \code{data.frame}, depends on whether 
#' the spatial grid is required.
#' 
#' @seealso \code{\link[UncerIn2]{Points-class}}, \code{\link[UncerIn2]{Grid.def}}, \code{\link[UncerIn2]{Grid.interpolation}}, \code{\link[UncerIn2]{uncertaintyInterpolation2-package}}
#' 
#' @name Grid.box
#' @docType methods
#' @rdname Grid.box
#' @aliases Grid.box,Points-method
#' @exportMethod Grid.box


setGeneric("Grid.box",
           function(object, ...)
           standardGeneric("Grid.box")
)


setMethod("Grid.box",
          signature(object = "Points"),
          definition = function(object, gridded = FALSE, cellsize = 10)
          {
            
            if (identical(gridded, FALSE)) 
            {
              x = seq(min(object@x), max(object@x),by = cellsize)
              y = seq(min(object@y), max(object@y),by = cellsize)
              
              grid = expand.grid (x,y)
              
              result=grid
              return(result)
            }
            
            else if (identical(gridded, TRUE))
            {
              x = seq(min(object@x), max(object@x),by = cellsize)
              y = seq(min(object@y), max(object@y),by = cellsize)
              
              grid = expand.grid (x,y)
              names(grid) = c("x", "y")
              gridded(grid) = ~x+y
              
              result=grid
              return(result) # class SpatialPixels 
            }   
            
            else {
              stop("argument `gridded` should be a logical TRUE or FALSE")
              
            }
    }
)


#' @title
#' Generation of grid from the results of interpolation functions
#'
#' @description
#' This method creates a grid from S4 object class \code{UncertainInterpolation}, defined in two possible 
#' ways - the user should decide if spatial grid is required (\code{gridded = TRUE/FALSE}). The size of 
#' the grid is recognized by the input data coordinates. User should define size of the single cell.
#'
#' @param object Input data type of S4 object class \code{UncertainInterpolation}.
#' @param gridded Logical value \code{TRUE/FALSE}. If true, then Spatial grid is made.
#' @param cellsize Numeric value for the size of cell.
#'
#' @usage
#' \S4method{Grid.interpolation}{UncertainInterpolation}(object,gridded,cellsize)
#' 
#' @return Returns an object of class \code{SpatialPixels} or \code{data.frame}, depends on whether 
#' the spatial grid is required.
#' 
#' @seealso \code{\link[UncerIn2]{UncertainInterpolation-class}}, \code{\link[UncerIn2]{Grid.box}}, \code{\link[UncerIn2]{Grid.def}}, \code{\link[UncerIn2]{uncertaintyInterpolation2-package}}
#' 
#' @name Grid.interpolation
#' @docType methods
#' @rdname Grid.interpolation
#' @aliases Grid.interpolation,UncertainInterpolation-method
#' @import sp
#' @exportMethod Grid.interpolation


setGeneric("Grid.interpolation",
           function(object, ...)
           standardGeneric("Grid.interpolation")
)


setMethod("Grid.interpolation",
          signature(object = "UncertainInterpolation"),
          definition = function(object, gridded = FALSE, cellsize = 10)
          {
            
            if (identical(gridded, FALSE)) 
            {
              x = seq(min(object@x), max(object@x),by = cellsize)
              y = seq(min(object@y), max(object@y),by = cellsize)
              
              grid = expand.grid (x,y)
              
              result=grid
              return(result)
            }
            
            else if (identical(gridded, TRUE))
            {
              x = seq(min(object@x), max(object@x),by = cellsize)
              y = seq(min(object@y), max(object@y),by = cellsize)
              
              grid = expand.grid (x,y)
              names(grid) = c("x", "y")
              gridded(grid) = ~x+y
              
              result=grid
              return(result) # class SpatialPixels 
            }     
            
            else {
              stop("argument `gridded` should be a logical TRUE or FALSE")
              
            }
    }
)