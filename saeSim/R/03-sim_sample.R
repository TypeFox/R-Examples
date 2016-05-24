#' Sampling component
#' 
#' One of the components which can be added to a \code{sim_setup}. This component can be used to add a sampling mechanism to the simulation set-up. A sample will be drawn after the population is generated (\code{\link{sim_gen}}) and variables on the population are computed (\code{\link{sim_comp_pop}}).
#' 
#' @param smplFun function which controls the sampling process.
#' @inheritParams sim_agg
#' 
#' @details Potentially you can define a \code{smplFun} yourself. Take care that it has one argument, named \code{dat} being the data as data.frame, and returns the sample as data.frame.
#' 
#' @seealso \code{\link{sample_number}}, \code{\link{sample_fraction}}
#' 
#' @export
#' @examples
#' # Simple random sample - 5% sample:
#' sim_base_lm() %>% sim_sample(sample_fraction(0.05))
#' 
#' # Simple random sampling proportional to size - 5% in each domain:
#' sim_base_lm() %>% sim_sample(sample_fraction(0.05, groupVars = "idD"))
#' 
#' # User defined sampling function:
#' sample_mySampleFun <- function(dat) {
#'   dat[sample.int(nrow(dat), 10), ]
#' }
#' 
#' sim_base_lm() %>% sim_sample(sample_mySampleFun)
sim_sample <- function(simSetup, smplFun = sample_number(size=5L, groupVars = "idD")) {
  sim_setup(simSetup, new("sim_fun", order = 5, call = match.call(), smplFun))
}