## This file is part of the UncertaintyInterpolation 2.0 package.
##
## Copyright 2015 Tomas Burian


#' @title
#' IDW interpolation
#'
#' @description
#' This function provides IDW interpolation over the input data enriched by the 
#' uncertainty. The input data must be an S4 object class of \code{UncertainPoints} 
#' and grid type of \code{Spatial}. Output object is type of S4 class 
#' \code{UncertainInterpolation}.
#'
#' @param object Input data. An object of \code{UncertainPoints} class.
#' @param grid Input Spatial grid.
#' @param nmax The number of nearest observations that should be used.
#' @param idp Inverse distance power.
#'
#' @usage
#' \S4method{idwUncertain}{UncertainPoints,Spatial}(object,grid,nmax,idp)
#' 
#' @return Returns an object of class \code{UncertainInterpolation}.
#' 
#' @seealso \code{\link[UncerIn2]{UncertainPoints-class}}, \code{\link[UncerIn2]{UncertainInterpolation-class}}, \code{\link[UncerIn2]{Grid.def}}, \code{\link[UncerIn2]{Grid.box}}, \code{\link[UncerIn2]{Grid.interpolation}}, \code{\link[gstat]{gstat}}, \code{\link[UncerIn2]{Plot}}, \code{\link[UncerIn2]{uncertaintyInterpolation2-package}}
#' 
#' @name idwUncertain
#' @docType methods
#' @rdname idwUncertain
#' @aliases idwUncertain,UncertainPoints,Spatial-method
#' @import gstat
#' @importFrom utils capture.output
#' @importFrom stats predict
#' @exportMethod idwUncertain
 

setGeneric("idwUncertain",
           function(object, grid, ...)
           standardGeneric("idwUncertain")
)


setMethod("idwUncertain",
          signature(object = "UncertainPoints", grid = "Spatial"),
          definition = function(object, grid, nmax = 7, idp = 2)           
          {
            
            if (!gridded(grid))
              return("Grid input data must be gridded!")
            
            data <- as.dataframe(object) # dataframe conversion
            
            a1 <- gstat(formula = uncertaintyLower ~ 1, locations = ~ x + y,
                        data = data, nmax = nmax, set = list(idp = idp))
            capture.output(a11 <- predict(a1, grid))
            
            a2 <- gstat(formula = modalValue ~ 1, locations = ~ x + y,
                        data = data, nmax = nmax, set = list(idp = idp))
            capture.output(a22 <- predict(a2, grid))
            
            a3 <- gstat(formula = uncertaintyUpper~ 1, locations = ~ x + y,
                        data = data, nmax = nmax, set = list(idp = idp))
            capture.output(a33 <- predict(a3, grid))
            
            new("UncertainInterpolation", x=grid$x, y=grid$y, uncertaintyLower=a11$var1.pred, 
                modalValue=a22$var1.pred, uncertaintyUpper=a33$var1.pred)
          }
)