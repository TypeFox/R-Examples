#' Death Rate Function
#'
#' Return the per capita rate of death (per unit time)
#'
#' @aliases deathFunc deathFunction
#'
#' @param stage (integer). The numbered life stage for which to return a death rate.
#' @param x Vector of state variables within the DDE solver. To access a variable use:
#' x$speciesName['stageName',strainNumber]
#' e.g. for species 'Bacteria', stage 'reproductive', strain 2 use
#' x$Bacteria['reproductive',2]
#' If there is only one stage and strain in species 'Food', for example, use
#' x$Food[1,1]
#' @param time (scalar). The current time point in the DDE solver.
#' @param species (integer). The numbered species for which to return a death rate.
#' @param strain (integer). The numbered strain for which to return a death rate.
#' @return per capita death rate (units of inverse time) for the strain, stage,species and time specified.
deathFuncDefault <- function(stage,x,time,species,strain){
  stop("You need to implement your own version of the death rate function.")
}
