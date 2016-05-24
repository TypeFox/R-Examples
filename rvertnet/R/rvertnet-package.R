#' Search VertNet archives using R
#'
#' There are a variety of ways to search VertNet
#' 
#' @section Search by term:
#' 
#' Search for _Aves_ in the state of _California_, limit to 10 records, e.g.:
#' 
#' \code{searchbyterm(class = "Aves", state = "California", lim = 10, verbose = FALSE)}
#' 
#' Search for _Mustela nigripes_ in the states of _Wyoming_ or _South Dakota_, 
#' limit to 20 records, e.g.:
#' 
#' \code{searchbyterm(genus = "Mustela", specificepithet = "nigripes", 
#'    state = "(wyoming OR south dakota)", limit = 20, verbose=FALSE)}
#' 
#' @section Big data:
#' Specifies a termwise search (like `searchbyterm()`), but requests that all available records 
#' be made available for download as a tab-delimited text file.
#' 
#' \code{bigsearch(genus = "ochotona", rf = "pikaRecords", email = "big@@search.luv")}
#' 
#' @section Spatial search:
#' \code{spatialsearch(lat = 33.529, lon = -105.694, radius = 2000, limit = 10, verbose = FALSE)}
#' 
#' @section Full text search:
#' Find records using a global full-text search of VertNet archives.
#' 
#' \code{vertsearch(taxon = "aves", state = "california")}
#' 
#' @importFrom methods is
#' @importFrom stats complete.cases setNames
#' @importFrom utils read.table
#' @importFrom httr GET content stop_for_status
#' @importFrom jsonlite fromJSON toJSON
#' @importFrom plyr compact rename
#' @importFrom dplyr rbind_all tbl_df src_sqlite tbl
#' @importFrom ggplot2 ggplot position_jitter aes geom_polygon 
#' geom_point labs theme_bw map_data
#' @import maps
#' @name rvertnet-package
#' @aliases rvertnet
#' @docType package
#' @keywords package
NULL
