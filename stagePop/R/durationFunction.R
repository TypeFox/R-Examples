#' Life stage duration function
#'
#' Return the duration of the life stages (time units). 
#'
#' @aliases durationFunc durationFunction
#'
#' @param stage (integer). The numbered life stage for which to return a duration
#' @param x Vector of state variables within the DDE solver. To access a particular variable use:
#' x$speciesName['stageName',strainNumber]
#' e.g. for species 'Bacteria', stage 'reproductive', strain 2 use
#' x$Bacteria['reproductive',2]
#' If there is only one stage and strain in species 'Food', for example, use
#' x$Food[1,1]
#' @param time (scalar). The current time point in the DDE solver.
#' @param species (integer). The numbered species for which to return a life stage duration.  
#' @param strain (integer). The numbered strain for which to return a life stage duration.  
#' @return Duration of the life stage for the stage, species, strain and time specified.
durationFuncDefault <- function(stage,x,time,species,strain){
  stop("You need to implement your own version of the duration function.")
}
