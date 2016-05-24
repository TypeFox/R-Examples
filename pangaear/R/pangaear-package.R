#' Client for the Pangaea Database
#'
#' @importFrom oai list_identifiers get_records list_sets list_metadataformats
#' list_identifiers id
#' @importFrom httr GET content stop_for_status write_disk config
#' @importFrom xml2 read_html xml_find_all xml_attr xml_text xml_find_one xml_parent
#' @importFrom tibble as_data_frame
#' @importFrom stats setNames
#' @importFrom utils head read.csv
#' @importFrom rappdirs user_cache_dir
#' @name pangaear-package
#' @aliases pangaear
#' @docType package
NULL
