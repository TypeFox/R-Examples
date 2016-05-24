#' Interface to the Dryad Solr API, and OAI-PMH service
#'
#' @importFrom httr GET content stop_for_status
#' @importFrom xml2 read_xml xml_find_all xml_ns xml_attr
#' @importFrom solr solr_search solr_facet solr_group solr_highlight
#' solr_mlt solr_stats
#' @importFrom oai id list_identifiers list_records list_metadataformats
#' list_sets get_records
#' @importFrom utils download.file
#' @name rdryad-package
#' @aliases rdryad
#' @docType package
NULL
