#' Package flora
#' 
#' Collect data from the Brazilian Flora checklist
#' (http://floradobrasil.jbrj.gov.br).
#' 
#' This package contains a set of tools solving problems that arise when one has
#' to collect taxonomic and distribution information for large datasets of
#' plants. Interacting with the Brazilian Flora Checklist website from a web browser is
#' often a slow and somewhat cumbersome process, especially if you are not sure
#' about the correct spelling of a name. With flora, however, you can:
#' \itemize{
#'   \item{get a suggestion for the correct spelling of a name from an incorrect one}
#'   \item{search for its current taxonomic status}
#'   \item{get its author(s), synonym(s), family, distribution, and lower taxa}
#'   \item{process lists of names and automatically solve synonyms and typing errors}
#' }
#' flora now holds all the data it needs. All functions can be used whilst offline.
#' @docType package
#' @import shiny httr dplyr
#' @importFrom "utils" "adist"
#' @name flora
#' @aliases flora flora-package
NULL