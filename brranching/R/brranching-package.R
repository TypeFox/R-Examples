#' Phylogenies from many sources
#'
#' @importFrom httr GET POST content stop_for_status
#' @importFrom curl curl_escape
#' @importFrom ape read.tree
#' @importFrom taxize tax_name
#' @importFrom phytools read.newick
#' @name brranching-package
#' @aliases brranching
#' @docType package
#' @author Scott Chamberlain <myrmecocystus@@gmail.com>
#' @keywords package
NULL

#' Lookup-table for family, genus, and species names for ThePlantList gymnosperms
#'
#' These names are from \url{http://www.theplantlist.org/}, collected
#' on 2015-11-11, and are from version 1.1 of their data. This data is
#' used in the function \code{\link{phylomatic_names}}.
#'
#' @format A data frame with 23,801 rows and 2 variables:
#' \describe{
#'   \item{family}{family name}
#'   \item{genus}{genus name}
#' }
#' @source \url{http://www.theplantlist.org/}
#' @name tpl
#' @docType data
#' @keywords data
NULL
