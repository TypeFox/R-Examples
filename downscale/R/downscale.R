################################################################################
# 
# downscale.R
# Version 1.4
# 08/05/2015
#
# Updates:
#   08/05/2015: extent now required
#   03/02/2015: output defined as class 'downscale'
#   03/02/2015: observed data included in output
#   02/02/2015: help file updated
#   02/02/2015: error check for model name and extent
#   30/01/2015: Thomas model added 
#
# Model area of occupancy against grain size for downscaling. 
#
# Args:
#   occupancies: data frame of observed area of occupancies and cell areas, or 
#                object from upgrain
#   model: function to use
#   extent: total area (same units as area)
#   tolerance: tolerance for integration of Thomas model. Lower numbers allow 
#              for greater accuracy but require longer processing times
#   starting_params: optional list of parameter values
#
# Returns:
#   list of three objects of class 'downscale'
#     $model    Downscaling model used
#     $pars     List of parameters estimated from optimisation procedure
#     $observed Dataframe of grain sizes and observed occupancies used for 
#               modelling
#
################################################################################

downscale <- function(occupancies,
                      model,
                      extent = NULL,
                      tolerance = 1e-6,
                      starting_params = NULL) {
  if(class(occupancies) == "upgrain") {
    extent <- occupancies$extent.stand
    occupancies <- occupancies$occupancy.stand[, -2]
  }
  
  if(class(occupancies) != "upgrain") {
    # error checking - input data frame correct
    if(ncol(occupancies) != 2) {
      stop("Input data must be a data frame with two columns (cell area and 
           occupancy")
    }
    
    # error checking - extent required
    if(is.null(extent)) {
      stop("Total extent required")
    }
    
    # error checking - extent larger than largest grain size
    if(extent < max(occupancies[, 2])) {
      stop("Total extent is smaller than the largest grain size! Are the units correct?")
    }
    
    # error checking - occupancies are between 0 and 1
    if(min(occupancies[, 2]) < 0) {
      stop("Occupancies must be proportion of cells occupied (values must be
           between 0 - 1)")
    }
    if(max(occupancies[, 2]) > 1) {
      stop("Occupancies must be proportion of cells occupied (values must be
           between 0 - 1)")
    }
  }
  
  # error checking - model name is correct
  if (model %in% c("Nachman", "PL", "Logis", "Poisson", "NB", "GNB", "INB",
                   "FNB", "Thomas") == FALSE) {
    stop("Model name invalid", call. = FALSE)
  }
  
  input.data <- DataInput(occupancy = occupancies[, 2],
                          area = occupancies[, 1],
                          extent = extent)
  model <- model
  if(is.null(starting_params)) {
    starting_params <- NULL
  }
  
  if ((model == "Nachman") | (model == "PL") | (model == "Logis") | 
      (model == "Poisson") | (model == "NB") | (model == "GNB") | 
      (model == "INB")){
    optim.pars <- suppressWarnings(OptimiseParameters(area =
                                                        input.data[!is.na(input.data[, "Occ"]),
                                                                   "Cell.area"], 
                                                      observed = 
                                                        input.data[!is.na(input.data[, "Occ"]),
                                                                   "Occ"],
                                                      model = model,
                                                      starting.params = starting_params))
  }
  
  if (model == "Logis") { 
    optim.pars <- suppressWarnings(
      OptimiseParametersLogis(area =
                                input.data[!is.na(input.data[,"Occ"]),
                                           "Cell.area"], 
                              observed = 
                                input.data[!is.na(input.data[,"Occ"]),
                                           "Occ"],
                              model = model,
                              starting.params = starting_params))
  }
  
  if (model == "GNB") { 
    optim.pars <- suppressWarnings(
      OptimiseParametersGNB(area =
                              input.data[!is.na(input.data[,"Occ"]),
                                         "Cell.area"], 
                            observed = 
                              input.data[!is.na(input.data[,"Occ"]),
                                         "Occ"],
                            model = model,
                            starting.params = starting_params))
  }
  
  if (model == "FNB") { 
    optim.pars <- suppressWarnings(
      OptimiseParametersFNB(area =
                              input.data[!is.na(input.data[,"Occ"]),
                                         "Cell.area"], 
                            observed = 
                              input.data[!is.na(input.data[,"Occ"]),
                                         "Occ"],
                            extent = extent,
                            model = model,
                            starting.params = starting_params))
  }
  
  if (model == "Thomas") { 
    optim.pars <- suppressWarnings(
      OptimiseParametersThomas(area =
                                 input.data[!is.na(input.data[,"Occ"]),
                                            "Cell.area"], 
                               observed = 
                                 input.data[!is.na(input.data[,"Occ"]),
                                            "Occ"],
                               extent = extent,
                               model = model,
                               tolerance = tolerance,
                               starting.params = starting_params))
  }
  observed <- data.frame("Cell.area" = input.data[,"Cell.area"],
                         "Occupancy" = input.data[,"Occ"])
  output <- list("model" = model,
                 "pars" = unlist(optim.pars),
                 "observed" = observed,
                 "extent" = extent)
  class(output) <- "downscale"
  return(output)
}
