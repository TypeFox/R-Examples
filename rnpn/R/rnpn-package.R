#' Interface to the National Phenology Network API
#'
#' @importFrom stats setNames
#' @importFrom data.table rbindlist setDF
#' @importFrom httr GET stop_for_status content
#' @importFrom jsonlite fromJSON
#' @importFrom plyr llply ldply ddply summarise rbind.fill
#' @name rnpn-package
#' @aliases rnpn
#' @docType package
#' @keywords package
NULL

#' Lookup-table for IDs of species and common names
#'
#' @name taxonlist
#' @docType data
#' @keywords data
#' @format A data.frame with 897 rows and 6 columns
#' \describe{
#'  \item{species_id}{species identifiers}
#'  \item{common_name}{common (vernacular) name}
#'  \item{genus}{genus name}
#'  \item{epithet}{epithet name}
#'  \item{itis_tsn}{ITIS taxonomic serial number (tsn)}
#'  \item{genus_epithet}{genus name + epithet name}
#' }
NULL
