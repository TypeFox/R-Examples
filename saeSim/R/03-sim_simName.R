#' Add a name to a sim_setup
#' 
#' Use this function to add a name to a \code{sim_setup} in case you are simulating different scenarios. This name will be added if you use the function \link{sim} for simulation
#' 
#' @inheritParams sim_agg
#' @param name a character
#' 
#' @export
#' @examples
#' sim_base_lm() %>% sim_simName("newName")
sim_simName <- function(simSetup, name) {
  slot(simSetup, "simName") <- name
  simSetup
}