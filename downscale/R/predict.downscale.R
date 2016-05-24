################################################################################
# 
# predict.downscale.R
# Version 1.2
# 30/01/2015
#
# Updates:
#   13/03/2015: if 0's predicted don't plot them
#   03/02/2015: plot function added
#   03/02/2015: output defined as class 'downscale'
#   03/02/2015: observed data included in output
#   02/02/2015: Help file added to
#   30/01/2015: Thomas model added
#
# Predict area of occupancy at fine grain sizes using parameters in a downscale
# object estimated from coarse grain sizes using downscale
#
# Args:
#   mod.fit: model output of class 'downscale' (created from function downscale)
#   new.areas: vector of grain sizes for model prediction
#   extent: total area (same units as newdata)- required only for FNB and Thomas
#           models
#   tolerance: tolerance for integration of Thomas model. Lower numbers allow
#              for greater accuracy but require longer processing times
#   plot: if TRUE plots observed and predicted occupancies against grain size on
#         a log-log plot
# 
# Returns:
#   list of three objects of class 'predict.downscale'
#     $model      Downscaling model used
#     $predicted  Dataframe of grain sizes and predicted occupancies
#     $observed   Dataframe of grain sizes and observed occupancies used for 
#                 modelling
#
################################################################################

predict.downscale <- function(object, 
                              new.areas, 
                              tolerance = 1e-6, 
                              plot = TRUE,
                              ...) {
  mod.fit <- object
  # error checking
  if (class(mod.fit) != "downscale"){
    stop("Input data not of class 'downscale'")
  }
  params <- as.list(mod.fit$pars)
  model <- mod.fit$model
  extent <- mod.fit$extent
  predict.function <- getFunction(paste("Predict", model, sep = ""))
  
  if ((model == "Nachman") | (model == "PL") | (model == "Logis") | 
        (model == "Poisson") | (model == "NB") | (model == "GNB") | 
        (model == "INB")){    
    AOO <- exp(predict.function(par = params, area = new.areas))
  }
  
  if (model == "FNB") {
    AOO <- exp(predict.function(par = params, 
                                area = new.areas, 
                                extent = extent))
  }
  
  if (model == "Thomas") {
    AOO <- exp(predict.function(par = params,
                                tolerance = tolerance,
                                area = new.areas, 
                                extent = extent))
  }
  expected <- data.frame("Cell.area" = new.areas,
                         "Occupancy" = AOO,
                         "AOO" = AOO * extent)
  output <- list("model" = model,
                 "predicted" = expected,
                 "observed" = mod.fit$observed)
  class(output) <- "predict.downscale"
  
  if (plot == TRUE) {
    par.original <- par()
    par.original <- list(mfrow = par.original$mfrow, mar = par.original$mar)
    par(mfrow = c(1, 1), mar = c(5, 5, 3, 1))
    
    plot.predict.downscale(output)
    
    par(mfrow = par.original$mfrow, mar = par.original$mar)
  }
  
  ### error checking in results
  if(sum(apply(expected, 1, function(x) is.na(x))[2 ,]) > 0) {
    warning("Predicted results may be innaccurate: one or more NA's predicted.")
  }
  
  if(sum(expected[, "Occupancy"] == 0, na.rm = TRUE) > 0) {
    warning("Predicted results may be innaccurate: one or more 0's predicted.")
  }
  
  if(sum(apply(expected, 1, function(x) is.na(x))[2 ,], na.rm = TRUE) > 0) {
    for(i in 1:(length(AOO) - 1)) {
      if(AOO[i] > AOO[i + 1]) {
        warning("Scaling is inconsistent: larger occupancies predicted at finer grain sizes.
               \nIf model = Thomas try a smaller tolerance value (e.g. 1e-7)")
      }
    }
  }
  
  return(output)
}