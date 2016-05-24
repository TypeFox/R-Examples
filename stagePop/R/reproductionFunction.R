#' Reproduction Function
#'
#' Return the rate of reproduction (amount per unit time)
#'
#' Unlike the other \code{\link{RateFunctions}} this rate function has
#' no stage argument as by definition it only pertains to the first
#' stage of life
#'
#' @aliases reproFunc reproductionFunction
#'
#' @param x Vector of state variables within the DDE solver. To access a variable use:
#' x$speciesName['stageName',strainNumber]
#' e.g. for species 'Bacteria', stage 'reproductive', strain 2 use
#' x$Bacteria['reproductive',2]
#' If there is only one stage and strain in species 'Food', for example, use
#' x$Food[1,1]
#' @param time (scalar). The current time point in the DDE solver.
#' @param species (integer). The numbered species for which to return a reproductive rate.
#' @param strain (integer). The numbered strain for which to return a rate.
#' @return Reproduction rate (amount per unit time) for the strain, species and time specified.
reproFuncDefault <- function(x,time,species,strain){
  stop("You need to implement your own version of the reproduction function.")
}
