#' Database of thermosensitive period of development for sex determination
#' @title Database of thermosensitive period of development for sex determination
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @docType data
#' @name TSP.list
#' @description Database of thermosensitive period of development for sex determination.\cr
#' This database can be used with the functions plot() or info.nests().\cr
#' The attributes TSP.begin.stages and TSP.end.stages for each dataframe give 
#' respectively the first and the last stages for TSP. Then the metrics for the limits of TSP 
#' are:\cr
#' TSP.list[[1]]$metric[TSP.list[[1]]$stages==attributes(TSP.list[[1]])$TSP.begin.stages] and \cr
#' TSP.list[[1]]$metric[TSP.list[[1]]$stages==1+attributes(TSP.list[[1]])$TSP.end.stages]
#' @references Mrosovsky, N., Pieau, C., 1991. Transitional range of temperature, pivotal 
#' temperatures and thermosensitive stages for sex determination in reptiles. 
#' Amphibia-Reptilia 12, 169-179.
#' @keywords datasets
#' @family Functions for temperature-dependent sex determination
#' @usage TSP.list
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(TSP.list)
#' names(TSP.list)
#' TSP.list[["Emys_orbicularis.mass"]]
#' attributes(TSP.list[["Emys_orbicularis.mass"]])$TSP.begin.stages
#' attributes(TSP.list[["Emys_orbicularis.mass"]])$TSP.end.stages
#' }
#' @format A list with dataframes including attributes
NULL
