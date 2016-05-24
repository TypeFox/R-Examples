## This file is part of the UncertaintyInterpolation 2.0 package.
##
## Copyright 2015 Tomas Burian


#' @title
#' Spline interpolation
#'
#' @description
#' This function provides Spline interpolation over the input data enriched by 
#' the uncertainty. The input data must be an S4 object class of 
#' \code{UncertainPoints} and grid type of \code{data.frame}. 
#' Output object is type of S4 class \code{UncertainInterpolation}.
#'
#' @param object Input data. An object of \code{UncertainPoints} class.
#' @param grid Input grid type of \code{dataframe}.
#' @param m A polynomial function of degree (m-1) will be included in 
#' the model as the drift (or spatial trend) component. Default is 
#' the value such that 2m-d is greater than zero where d is the dimension of x.
#' @param p Polynomial power for Wendland radial basis functions. Default is 
#' 2m-d where d is the dimension of x.
#' @param scale.type The independent variables and knots are scaled to the specified 
#' scale.type. By default the scale type is "range", whereby the locations are 
#' transformed to the interval (0,1) by forming (x-min(x))/range(x) for each x. 
#' Scale type of "user" allows specification of an x.center and x.scale by the user. 
#' The default for "user" is mean 0 and standard deviation 1. Scale type of "unscaled" 
#' does not scale the data.
#' @param lon.lat If TRUE locations are interpreted as lognitude and latitude and great circle distance is used to find distances among locations.
#' @param miles If TRUE great circle distances are in miles if FALSE distances are in kilometers.
#' @param method Determines what "smoothing" parameter should be used. The default is to estimate standard GCV Other choices are: GCV.model, GCV.one, RMSE, pure error and REML. The differences are explained in the Krig help file.
#' @param GCV If TRUE the decompositions are done to efficiently evaluate the estimate, GCV function and likelihood at multiple values of lambda.
#'
#' @usage
#' \S4method{splineUncertain}{UncertainPoints,data.frame}(object, grid, m = NULL, p = NULL, 
#'    scale.type = "range", lon.lat = FALSE, miles = TRUE, method = "GCV", GCV = TRUE)
#' 
#' @return Returns an object of class \code{UncertainInterpolation}.
#' 
#' @seealso \code{\link[UncerIn2]{UncertainPoints-class}}, \code{\link[UncerIn2]{UncertainInterpolation-class}}, \code{\link[UncerIn2]{Grid.def}},\code{\link[UncerIn2]{Grid.box}}, \code{\link[UncerIn2]{Grid.interpolation}}, \code{\link[fields]{Tps}}, \code{\link[UncerIn2]{Plot}}, \code{\link[UncerIn2]{uncertaintyInterpolation2-package}}
#' 
#' @name splineUncertain
#' @docType methods
#' @rdname splineUncertain
#' @aliases splineUncertain,UncertainPoints,data.frame-method
#' @import fields
#' @importFrom stats predict
#' @exportMethod splineUncertain


setGeneric("splineUncertain",
           function(object, grid, ...)
           standardGeneric("splineUncertain")
          )


setMethod("splineUncertain",
          signature(object = "UncertainPoints", grid = "data.frame"),
          definition = function(object, grid, m = NULL, p = NULL, scale.type = "range", lon.lat = FALSE,
                                miles = TRUE, method = "GCV", GCV = TRUE)
          {
            
            if(!(inherits(grid, "data.frame"))){
              stop( paste("Grid: ", deparse(substitute(testData))," is not of type data.frame." , sep=""))
            }
            
            a = as.dataframe(object)
            loc = cbind(a$x, a$y)
            
            a1 <- Tps(loc, object@uncertaintyLower, m = m, p = p, scale.type = scale.type, 
                      lon.lat = lon.lat, miles = miles, method = method, GCV = GCV)
            a11 <- predict(a1, grid)

            a2 <- Tps(loc, object@modalValue, m = m, p = p, scale.type = scale.type, 
                      lon.lat = lon.lat, miles = miles, method = method, GCV = GCV) 
            a22 <- predict(a2, grid)

            a3 <- Tps(loc, object@uncertaintyUpper, m = m, p = p, scale.type = scale.type, 
                      lon.lat = lon.lat, miles = miles, method = method, GCV = GCV) 
            a33 <- predict(a3, grid)
            
            new("UncertainInterpolation", x=grid$x, y=grid$y, uncertaintyLower=a11[,1], modalValue=a22[,1], uncertaintyUpper=a33[,1])
          }
)