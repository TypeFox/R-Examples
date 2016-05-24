#' DataCite R client.
#'
#' @section OAI-PMH functions:
#'
#' \itemize{
#'   \item \code{\link{dc_oai_getrecord}} - Get records
#'   \item \code{\link{dc_oai_identify}} - identify the OAI-PMH service
#'   \item \code{\link{dc_oai_listidentifiers}} - List identifiers
#'   \item \code{\link{dc_oai_listmetadataformats}} - List metadata formats
#'   \item \code{\link{dc_oai_listrecords}} - List records
#'   \item \code{\link{dc_oai_listsets}} - List sets
#' }
#'
#' @section Search functions:
#'
#' \itemize{
#'   \item \code{\link{dc_search}} - General search, only returns documents
#'   \item \code{\link{dc_facet}} - Faceting only (w/o general search)
#'   \item \code{\link{dc_mlt}} - More like this (w/o general search)
#'   \item \code{\link{dc_stats}} - Stats search (w/o general search)
#' }
#'
#' @section Vignettes:
#'
#' Coming soon...
#'
#' @importFrom solrium solr_search solr_facet solr_stats solr_mlt solr_connect solr_settings
#' @importFrom oai get_records id list_metadataformats list_records list_sets
#' @name rdatacite-package
#' @aliases rdatacite
#' @docType package
#' @author Scott Chamberlain \email{myrmecocystus@@gmail.com}
#' @keywords package
NULL
