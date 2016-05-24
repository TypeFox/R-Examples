#' Reproduce a result
#' @description Helps reproduce the result of a simulation by restoring the RNG to the state of a supplied null model object.
#' @param model the model object containing the result to be reproduced.
#' @details Works by resetting the RNG state to what it was for a given EcoSimR simulation.  This only works if the user saved the seed with the saveSeed parameter.
#' @examples \dontrun{
#' finchMod <- cooc_null_model(dataWiFinches, algo="sim1",saveSeed=T)
#' ## Check model output
#' mean(finchMod$Sim)
#' 
#' reproduce_model(finchMod)
#' 
#' finchMod <- cooc_null_model(dataWiFinches, algo="sim1")
#' ## Check model output is the same as before
#' mean(finchMod$Sim)
#' reproduce_model(finchMod)
#' }
#' @export

reproduce_model <- function(model) {
  if (model$Reproducible){
    assign(".Random.seed", model$RandomSeed, .GlobalEnv)
  } else {
    stop("You didn't save the seed for this model run, please run again and set saveSeed = TRUE")
  }
  
}
