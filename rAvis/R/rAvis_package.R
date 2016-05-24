#' We developed several functions to explore and donwload 
#' the information stored in ProyectoAVIS database (www.proyectoavis.com), 
#' in an easy and visual way.
#' 
#' We programmed two main functions to set flexible queries 
#' about the species occurrences and the birdwatcher 
#' observations: avisQuerySpecies and avisQueryContributor. 
#' Besides, there are also general functions 
#' to explore the database, like avisMapSpecies.
#' 
#' @name rAvis
#' @aliases rAvis-package
#' @docType package
#' @title rAvis: An R-package to download the information stored in Proyecto AVIS, 
#' a citizen science bird project.
#' @author Javier Gonzalez Hernandez \email{javigzz@@yahoo.es}
#' @author Sara Varela \email{svarela@@paleobiogeography.org}
#' @references Varela S, Gonzalez-Hernandez J, Casabella E, Barrientos R (2014) 
#' rAvis: An R-Package for Downloading Information Stored in Proyecto AVIS, 
#' a Citizen Science Bird Project. PLoS ONE 9(3): e91650. doi: 10.1371/journal.pone.0091650
#' @keywords package
#' @details \tabular{ll}{
#' Package: \tab rAvis \cr
#' Type: \tab Package\cr
#' Version: \tab 0.1.3 \cr
#' Date: \tab 2015-06-17\cr
#' License: \tab GPL-2 \cr
#' }
#' 
#'@seealso {
#' http://proyectoavis.com
#' }
#'@examples \dontrun{
#' avisSpeciesSummary()
#'
#' avisMapSpecies ("Pica pica", maptype="phys")
#'
#' avisQuerySpecies(list("Bubo bubo", "Tyto alba"), args = list(year = 2012))
#' }
#'
NULL