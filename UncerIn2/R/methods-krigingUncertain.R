## This file is part of the UncertaintyInterpolation 2.0 package.
##
## Copyright 2015 Tomas Burian


#' @title
#' Krigging and Fuzzy kriging interpolation
#'
#' @description
#' This function provides kriging interpolation over the input data enriched by the 
#' uncertainty with multiple dispatch. The input data must be an S4 object class 
#' of \code{UncertainPoints} or \code{Points} in case 
#' of Fuzzy kriging. Grid type of \code{Spatial}. Output object is type of S4 class 
#' \code{UncertainInterpolation} or in case of Fuzzy kriging S4 class  
#' \code{FuzzyInterpolation}. For more informations about fuzzy kriging 
#' interpolation read Details.
#'
#' @param object Input data. An object of \code{UncertainPoints} class.
#' @param grid Input Spatial grid.
#' @param ... Additional arguments to be passed to \code{f}.
#' @param data_variogram An optional way to provide a different dataset for the building of the variogram then for the spatial interpolation.
#' @param block Use this parameter to pass on a specification for the block size. e.g. c(1000,1000)
#' @param model List of models that will be tested during automatic variogram fitting
#' @param kappa List of values for the smoothing parameter of the Matern model that will be tested during automatic variogram fitting.
#' @param fix.values Can be used to fix a variogram parameter to a certain value. It consists of a list with a length of three. The items describe the fixed value for the nugget, range and sill respectively. Setting the value to NA means that the value is not fixed. Is passed on to autofitVariogram.
#' @param remove_duplicates logical, remove duplicates from input
#' @param verbose logical, if TRUE autoKrige will give extra information on the fitting process
#' @param GLS.model If a variogram model is passed on through this parameter a Generalized Least Squares sample variogram is calculated.
#' @param start_vals Can be used to give the starting values for the variogram fitting. The items describe the fixed value for the nugget, range and sill respectively.
#' @param miscFitOptions Additional options to set the behavior of autofitVariogram.
#' 
#' @param krigingModel Model of the kriging used in the calculations.
#' @param psills Defined psills.
#' @param ranges Defined ranges.
#' @param nuggets Defined nuggets.
#' @param models Variogram model of dependent variable (or its residuals).
#' @param vgm_start Modal variogram selected for the data.
#' @param logResults Was the dataset logaritmized prior to the calculation?
#'
#' @usage
#' \S4method{krigingUncertain}{UncertainPoints,Spatial}(object, grid, data_variogram = data, block = 0, 
#'    model = c("Sph", "Exp", "Gau", "Ste"), kappa = c(0.05, seq(0.2, 2, 0.1), 5, 10), 
#'    fix.values = c(NA,NA,NA), remove_duplicates = TRUE, verbose = FALSE, GLS.model = NA, 
#'    start_vals = c(NA,NA,NA), miscFitOptions = list())
#' 
#' \S4method{krigingUncertain}{Points,Spatial}(object, grid, krigingModel, psills, ranges, nuggets, models, 
#'    vgm_start, logResults=FALSE)
#' 
#' @return Returns an object of class \code{UncertainInterpolation} or \code{FuzzyInterpolation}.
#' 
#' @seealso \code{\link[UncerIn2]{UncertainPoints-class}}, \code{\link[UncerIn2]{Points-class}}, \code{\link[UncerIn2]{UncertainInterpolation-class}}, \code{\link[UncerIn2]{FuzzyInterpolation-class}}, \code{\link[UncerIn2]{Grid.def}},\code{\link[UncerIn2]{Grid.box}}, \code{\link[UncerIn2]{Grid.interpolation}}, \code{\link[automap]{autoKrige}}, \code{\link[UncerIn2]{Plot}}, \code{\link[UncerIn2]{uncertaintyInterpolation2-package}}, \url{https://github.com/JanCaha/Hais2015-paper}
#' 
#' @details The function for Fuzzy kriging and its processes were taken from source code Jan Caha. For more informations and details visit \url{https://github.com/JanCaha/Hais2015-paper}.
#' 
#' @name krigingUncertain
#' @docType methods
#' @rdname krigingUncertain
#' @aliases krigingUncertain,UncertainPoints,Spatial-method
#'          krigingUncertain,Points,Spatial-method
#' @import automap
#' @exportMethod krigingUncertain
 

setGeneric("krigingUncertain",
           function(object, grid, ...)
           standardGeneric("krigingUncertain")
)


setMethod("krigingUncertain",
          signature(object = "UncertainPoints", grid = "Spatial"),
          definition = function(object, grid, data_variogram = data, block = 0, model = c("Sph", "Exp", "Gau", "Ste"), 
                                kappa = c(0.05, seq(0.2, 2, 0.1), 5, 10), fix.values = c(NA,NA,NA), remove_duplicates = TRUE, 
                                verbose = FALSE, GLS.model = NA, start_vals = c(NA,NA,NA), miscFitOptions = list())           
          {
            
            if (!gridded(grid))
            return("Grid input data must be gridded!")
            
            data <- as.dataframe(object) # dataframe conversion
            coordinates(data) <- ~ x+y
            
            capture.output(a1 <- autoKrige(uncertaintyLower~x+y, data, grid, data_variogram = data_variogram, block = block, 
                                           model = model, kappa = kappa, fix.values = fix.values, remove_duplicates = remove_duplicates,
                                           verbose = verbose, GLS.model = GLS.model, start_vals = start_vals, miscFitOptions = miscFitOptions))
            a11 = a1$krige_output
            
            capture.output(a2 <- autoKrige(modalValue~x+y, data, grid, data_variogram = data_variogram, block = block, 
                                          model = model, kappa = kappa, fix.values = fix.values, remove_duplicates = remove_duplicates,
                                          verbose = verbose, GLS.model = GLS.model, start_vals = start_vals, miscFitOptions = miscFitOptions))
            a22 = a2$krige_output
            
            capture.output(a3 <- autoKrige(uncertaintyUpper~x+y, data, grid, data_variogram = data_variogram, block = block, 
                                           model = model, kappa = kappa, fix.values = fix.values, remove_duplicates = remove_duplicates,
                                           verbose = verbose, GLS.model = GLS.model, start_vals = start_vals, miscFitOptions = miscFitOptions))
            a33 = a3$krige_output
        
            new("UncertainInterpolation", x=grid$x, y=grid$y, uncertaintyLower=a11$var1.pred, modalValue=a22$var1.pred, uncertaintyUpper=a33$var1.pred)
          }
)


setMethod("krigingUncertain",
          signature(object = "Points", grid = "Spatial"),
          definition = function(object, grid, krigingModel, psills, ranges, nuggets, models, 
                                vgm_start, logResults=FALSE)
          {
            
            if (!gridded(grid))
              return("Grid input data must be gridded!")
            
            data <- as.dataframe(object) # dataframe conversion
            coordinates(data) <- ~ x+y
            
            numberOfCalculations = length(psills) * length(ranges) * length(nuggets) * length(models)
            
            # initial prediction is done based on modal model
            pred = krige(krigingModel, data, grid, model = vgm_start)
            
            # were the data logaritmized prior to the calculation?
            if(logResults){
              min_data = exp(pred$var1.pred)
              max_data = exp(pred$var1.pred)
            }
            else{
              min_data = pred$var1.pred
              max_data = pred$var1.pred
            }
            
            calculationNumber = 1
            
            # for all combinations of the sills, ranges, nuggers and also models(if there is more than one)
            for (psill in psills){
              for (range in ranges){
                for(nugget in nuggets){
                  for(model in models){
                    
                    #prepare the varigogram
                    vgm <- vgm(psill = psill, model= model,range= range, nugget=nugget)
                    
                    #calculated the krigging
                    pred <- krige(krigingModel, data, grid, model = vgm_start)
                    
                    #if the dataset was logaritmized we need to exponentiate the outcome
                    if(logResults){
                      temp_data = exp(pred$var1.pred)
                    }
                    else{
                      temp_data = pred$var1.pred
                    }
                    
                    #compare the obtained preditions to minimal and maximal values that we allready have
                    #if the value is outside of the range, than adjust the range
                    for(i in 1:length(temp_data)){
                      if(temp_data[i] < min_data[i]){
                        min_data[i] = temp_data[i]
                      }
                      
                      if(max_data[i] < temp_data[i]){
                        max_data[i] = temp_data[i]
                      }
                    }
                    
                    print(paste("Calculation",calculationNumber,"of",numberOfCalculations,"done.",(calculationNumber/numberOfCalculations)*100,"%.", sep = " "))
                    calculationNumber = calculationNumber + 1
                  }
                } 
              }
            }
            
            predModal = krige(z~1, data, grid, model = vgm_start)
            b1 = predModal$var1.pred
            
            #prepare the result and return it from the function
            new("FuzzyInterpolation", x=grid$x, y=grid$y, minimal=min_data, modalValue=b1, maximal=max_data)
          }
)