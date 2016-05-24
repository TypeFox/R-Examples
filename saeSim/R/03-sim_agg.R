#' Aggregation component
#' 
#' One of the components which can be added to a simulation set-up. Aggregating the data is a simulation component which can be used to aggregate the population or sample. The aggregation will simply be done after the sampling, if you haven't specified any sampling component, the population is aggregated (makes sense if you draw samples directly from the model).
#' 
#' @param simSetup a \code{sim_setup}.
#' @param aggFun function which controls the aggregation process. At the moment only \code{\link{agg_all}} is defined. 
#' 
#' @details Potentially you can define an \code{aggFun} yourself. Take care that it only has one argument, named \code{dat}, and returns the aggregated data as \code{data.frame}.
#' 
#' @seealso \code{\link{agg_all}}, \code{\link{sim_gen}}, \code{\link{sim_comp_pop}}, \code{\link{sim_sample}}, , \code{\link{sim_comp_sample}}
#' 
#' @export
#' @examples
#' # Aggregating the population:
#' sim_base_lm() %>% sim_agg()
#' 
#' # Aggregating after sampling:
#' sim_base_lm() %>% sim_sample() %>% sim_agg()
#' 
#' # User aggFun:
#' sim_base_lm() %>% sim_agg(function(dat) dat[1, ])
sim_agg <- function(simSetup, aggFun = agg_all()) {
  sim_setup(simSetup, new("sim_fun", order = 7, call = match.call(), aggFun))
}