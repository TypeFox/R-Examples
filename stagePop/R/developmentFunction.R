#' Development Function
#'
#' Return the rate of development (per unit time)
#'
#' This function will only be called
#' when the \code{timeDependDuration} argument to
#' \code{\link{popModel}} contains TRUE values, otherwise the
#' development rate within a stage is irrelevant to the model. 
#'
#' @aliases develFunc developmentFunction
#'
#' @param stage (integer). The numbered life stage for which to return a developement rate.
#' @param x Vector of state variables within the DDE solver. To access a variable use:
#' x$speciesName['stageName',strainNumber]
#' e.g. for species 'Bacteria', stage 'reproductive', strain 2 use
#' x$Bacteria['reproductive',2]
#' If there is only one stage and strain in species 'Food', for example, use
#' x$Food[1,1]
#' @param time (scalar). The current time point in the DDE solver.
#' @param species (integer). The numbered species for which to return a development rate.
#' @param strain (integer). The numbered strain for which to return a developement rate.
#' @return Development rate (units of inverse time) for the strain, stage, species and time specified.
develFuncDefault <- function(stage,x,time,species,strain){
  stop("You need to implement your own version of the development function.")
}
  
