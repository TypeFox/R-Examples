################################################################################
# 
# ensemble.downscale.R
# Version 1.2
# 13/03/2015
#
# Updates:
#   13/03/2015: if 0's predicted don't plot them
#   24/02/2015 means calculated as mean of log occupancies
#   24/02/2015 warning messages for inconsistent results
#   23/02/2015 different tolerances allowed for modelling and predicting
#   05/02/2015 plot = TRUE argument added
#   05/02/2015 improved warning control
#
# Ensemble modelling
#
# Args:
#   occupancies: data frame of observed area of occupancies and cell areas, or 
#                object from upgrain
#   new.areas: vector of grain sizes (in same units as occupancy) for which area 
#               of occupancy will be predicted.
#   extent: total area in same units as occupancy
#   tolerance_mod: The tolerance
#           used during integration in the Thomas model during optimisation of
#           parameters.
#   tolerance_pred: The tolerance used during the prediction stage. 
#   tolerance_hui: The tolerance used in the Hui model
#   models: vector of chosen downscaling models. Default models = "all" runs all
#           available models.
#   plot: if TRUE predictions of all models are plotted against grain size along
#         with the mean of all models.
#   verbose: if TRUE prints updates on modelling status.
# Returns:
#    a dataframe. The first column cell.area is the grain sizes used for 
#    predictions. The final column Means are the mean predictions of all models
#    for each grain size. Intermediate columns are the predicted occupancies for
#    the selected downscaling models.
#
################################################################################

ensemble.downscale <- function(occupancies,
                               new.areas,
                               extent = NULL,
                               cell.width = NULL,
                               models = "all",
                               tolerance_mod = 1e-6,
                               tolerance_pred = 1e-6,
                               tolerance_hui = 1e-6,
                               starting_params = NULL,
                               plot = TRUE,
                               verbose = TRUE) {
  ## error checking - model name is correct
  apply(as.data.frame(models), 1, function(x)
    if (x %in% c("Nachman", "PL", "Logis", "Poisson", "NB", "GNB", "INB",
                 "FNB", "Thomas", "Hui", "all") == FALSE) {
      stop("Model name invalid", call. = FALSE)
    })
  
  if (length(models) == 1) {
    if (models == "all") {
      model.list <- c("Nachman","PL","Logis","Poisson","NB",
                      "GNB","INB","FNB","Thomas", "Hui")
    }
    if (models != "all") {
      stop("Only one model selected: ensemble modelling not applicable", 
           call. = FALSE)
    }
  }
  
  ## error checking - if not an upgrain object
  if(class(occupancies) != "upgrain") {
    # error checking - extent given
    if(is.null(extent)) {
      stop("No extent given (occupancies is not of class 'upgrain'")
    }
    
    # error checking - extent larger than largest grain size
    if(extent < max(occupancies[, 2])) {
      stop("Total extent is smaller than the largest grain size! Are the units correct?")
    }
    
    # error checking - for hui model occupancies is upgrain object
    if (length(models) == 1) {
      if(models == "all") {
       stop("Hui model can not be run if occupancies is not of class 'upgrain'")
      }
    }
    if(sum(models == "Hui") > 0) {
      stop("Hui model can not be run if occupancies is not of class 'upgrain'")
    }
  }
  
  if(class(occupancies) == "upgrain") {
    cell.width <- raster::res(occupancies$atlas.raster.stand)[1]
    extent <- occupancies$extent.stand
  }
  
  # error checking - if model = Hui then cell.width must be present
  if(sum(models == "Hui") > 0) {
    if(is.null(cell.width)) {
      stop("Cell.width must be specified for the Hui model")
    }
  }
  
  ## data handling
  if (length(models) > 1) {
    model.list <- models
  }
  
  starting_params_opts <- starting_params
  starting_params_mods <- NULL
  if(!is.null(starting_params)) {
    starting_params_mods <- names(starting_params_opts)
  }
  
  all.predicted <- as.data.frame(matrix(NA, 
                                        ncol = (length(model.list) + 1),
                                        nrow = length(new.areas)))
  colnames(all.predicted) <- c("Cell.area", model.list)
  all.predicted[, "Cell.area"] <- new.areas
  
  # modelling
  for (i in 1:length(model.list)) {
    model.run<-model.list[i]
    
    ## see if there are user-inputted starting parameters
    if(sum(starting_params_mods == model.run) == 1) {
      starting_params <- starting_params_opts[[model.run]]
    } else {
      starting_params <- NULL
    }
    
    if(verbose == TRUE){
      cat(paste(model.run, "model is running..."))
    }
    
    if(model.run != "Hui") {
      mod <- downscale(occupancies = occupancies,
                       model = model.run,
                       extent = extent,
                       tolerance = tolerance_mod,
                       starting_params = starting_params)
      est <- predict.downscale(object = mod,
                     new.areas = new.areas,
                     extent = extent,
                     tolerance = tolerance_pred,
                     plot = FALSE)
      all.predicted[, i + 1] <- est$predicted[, "Occupancy"]
    }
    if(model.run == "Hui") {
      new.areas.hui <- new.areas[new.areas < (cell.width ^ 2)]
      est <- hui.downscale(atlas.data = occupancies, 
                           cell.width = cell.width, 
                           new.areas = new.areas.hui, 
                           extent = extent,
                           tolerance = tolerance_hui)
      all.predicted[1:length(new.areas.hui), 
                    i + 1] <- est$predicted[, "Occupancy"]
    }
    
    if(verbose == TRUE){
      cat(paste("  complete", "\n"))
    }
  }
  all.predicted$Means <- exp(rowMeans(log(all.predicted[, -1]), na.rm = TRUE))
  aoo.predicted <- all.predicted
  aoo.predicted[, -1] <- aoo.predicted[, -1] * extent
  
  
  # plotting
  if (plot == TRUE) {
    par.original <- par()
    par.original <- list(mfrow = par.original$mfrow, mar = par.original$mar)
    par(mfrow = c(3, ceiling(length(model.list) / 3)), mar = c(5, 5, 3, 1))
    
    for (i in 1:length(model.list)) {
      predicted <- all.predicted
      predicted[predicted == 0] <- NA
      plot(predicted[, 2] ~ all.predicted[, "Cell.area"],
           type = "n",
           log = "xy",
           xlim = c(min(c(all.predicted[, "Cell.area"],
                          est$observed[, "Cell.area"]), na.rm = TRUE),
                    max(c(all.predicted[, "Cell.area"],
                          est$observed[, "Cell.area"]), na.rm = TRUE)),
           ylim = c(min(predicted[, -1], na.rm = TRUE), 1),
           xlab = "Log cell area",
           ylab = "Log occupancy",
           main = paste(model.list[i], "model"))
      
      points(est$observed[, "Occupancy"] ~ est$observed[, "Cell.area"],
             type="b",
             lwd=2)
      points(predicted[, "Means"] ~ all.predicted[, "Cell.area"],
             type="b",
             lwd=2,
             col = "dark grey")
      points(predicted[, i + 1] ~ all.predicted[, "Cell.area"],
             type="b",
             lwd=2,
             col = "red")
    }
    par(mfrow = par.original$mfrow, mar = par.original$mar)
  }
  output <- list(Occupancy = all.predicted,
                 AOO = aoo.predicted)
  return(output)
}

