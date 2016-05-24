## This file is part of the UncertaintyInterpolation 2.0 package.
##
## Copyright 2015 Tomas Burian


#' @title
#' Estimating of variograms
#'
#' @description
#' This function can estimate variogram over the input data of S4 object 
#' class \code{UncertainPoints}. It is a bonus function of this whole 
#' package - just for user comfort. Output object is type of S4 class \code{Variogram} 
#' which can be printed by function \code{show} from UncerIn2 package. 
#' It also visualizes automaticaly all estimated variograms during the process. 
#'
#' @param object Input data type of S4 object class UncertainPoints.
#' @param model The list of variogrammodels that will be tested.
#' @param kappa Smoothing parameter of the Matern model. Provide a list if you want to check more than one value.
#' @param fix.values Can be used to fix a variogram parameter to a certain value. It consists of a list with a length of three. The items describe the fixed value for the nugget, range and sill respectively.
#' @param verbose Logical, if TRUE the function will give extra feedback on the fitting process.
#' @param GLS.model If a variogram model is passed on through this parameter a Generalized Least Squares sample variogram is calculated.
#' @param start_vals Can be used to give the starting values for the variogram fitting. The items describe the fixed value for the nugget, range and sill respectively. They need to be given in that order. Setting the value to NA means that the value will be automatically chosen.
#' @param miscFitOptions A list with named arguments that provide additional control over the fitting process. For example: list(merge.small.bins = TRUE). If the list is empty, autofitVariogram uses default values.
#' 
#' @usage
#' \S4method{variogram}{UncertainPoints}(object, model = c("Sph", "Exp", "Gau", "Ste"), 
#'    kappa = c(0.05, seq(0.2, 2, 0.1), 5, 10), fix.values = c(NA,NA,NA), verbose = FALSE, 
#'    GLS.model = NA, start_vals = c(NA,NA,NA), miscFitOptions = list())
#'
#' @return Returns an object of class \code{Variogram}.
#'
#' @details For the estimation of variogram were used functions from package automap.
#' 
#' @seealso \code{\link[UncerIn2]{UncertainPoints-class}},  \code{\link[automap]{autofitVariogram}}, \code{\link[UncerIn2]{uncertaintyInterpolation2-package}}
#'
#' @name variogram
#' @docType methods
#' @rdname variogram
#' @aliases variogram,UncertainPoints-method
#' @exportMethod variogram


setGeneric("variogram",
           function(object, ...)
             standardGeneric("variogram")
)


setMethod("variogram",
          signature(object = "UncertainPoints"),
          definition = function(object, model = c("Sph", "Exp", "Gau", "Ste"), kappa = c(0.05, seq(0.2, 2, 0.1), 5, 10), 
                                fix.values = c(NA,NA,NA), verbose = FALSE, GLS.model = NA, start_vals = c(NA,NA,NA), 
                                miscFitOptions = list())         
          {
            
            data <- as.dataframe(object) # dataframe conversion
            coordinates(data) <- ~ x+y # coordinates
            
            a1 = autofitVariogram(uncertaintyLower~x+y, data, model = model, kappa = kappa, fix.values = fix.values, 
                                  verbose = verbose, GLS.model = GLS.model, start_vals = start_vals, miscFitOptions = miscFitOptions)
            
            a2 = autofitVariogram(modalValue~x+y, data, model = model, kappa = kappa, fix.values = fix.values, 
                                  verbose = verbose, GLS.model = GLS.model, start_vals = start_vals, miscFitOptions = miscFitOptions)
            
            a3 = autofitVariogram(uncertaintyUpper~x+y, data, model = model, kappa = kappa, fix.values = fix.values, 
                                  verbose = verbose, GLS.model = GLS.model, start_vals = start_vals, miscFitOptions = miscFitOptions)
          
            plot(a1)
            plot(a2)
            plot(a3)
            
            a11 = new("Vario", model = as.character(a1$var_model$model), psill = a1$var_model$psill, range = a1$var_model$range,
                kappa = a1$var_model$kappa, ang1 = a1$var_model$ang1, ang2 = a1$var_model$ang2, ang3 = a1$var_model$ang3,
                anis1 = a1$var_model$anis1, anis2 = a1$var_model$anis2)
            
            a22 = new("Vario", model = as.character(a2$var_model$model), psill = a2$var_model$psill, range = a2$var_model$range,
                kappa = a2$var_model$kappa, ang1 = a2$var_model$ang1, ang2 = a2$var_model$ang2, ang3 = a2$var_model$ang3,
                anis1 = a2$var_model$anis1, anis2 = a2$var_model$anis2)            
            
            a33 = new("Vario", model = as.character(a3$var_model$model), psill = a3$var_model$psill, range = a3$var_model$range,
                kappa = a3$var_model$kappa, ang1 = a3$var_model$ang1, ang2 = a3$var_model$ang2, ang3 = a3$var_model$ang3,
                anis1 = a3$var_model$anis1, anis2 = a3$var_model$anis2)
            
            new("Variogram", minimal = a11, modalValue = a22, maximal = a33)
          }
)

