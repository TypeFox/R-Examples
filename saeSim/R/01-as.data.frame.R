#' as.data.frame method for sim_setup
#' 
#' Use this method to get a single simulated data.frame out of a sim_setup object.
#' 
#' @param x a sim_setup
#' @param row.names will have no effect
#' @param optional will have no effect
#' @param ... will have no effect
#' 
#' @export
#' @method as.data.frame sim_setup
as.data.frame.sim_setup <- function(x, row.names = NULL, optional = FALSE, ...) {
  sim_run_once(x) %>% as.data.frame(row.names = row.names, optional = optional, ...)
}
