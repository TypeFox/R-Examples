#' Retrieve record headers
#'
#' @export
#'
#' @param from specifies that records returned must have been created/update/deleted
#'     on or after this date.
#' @param until specifies that records returned must have been created/update/deleted
#'     on or before this date.
#' @param set specifies the set that returned records must belong to.
#' @param prefix specifies the metadata format that the records will be
#'     returned in. One of: oai_dc (default), oai_datacite, or datacite. See
#'     \code{Prefixes} for more info.
#' @param token a token previously provided by the server to resume a request
#'     where it last left off.
#' @param ... Curl options passed on to \code{\link[httr]{GET}}
#' @inheritParams dc_oai_listmetadataformats
#' @examples \dontrun{
#' today <- format(Sys.Date(), "%Y-%m-%d")
#' dc_oai_listidentifiers(from = today)
#' dc_oai_listidentifiers(from = '2011-06-01T', until = '2011-07-01T')
#' dc_oai_listidentifiers(set = "REFQUALITY")
#' }
dc_oai_listidentifiers <- function(from = NULL, until = NULL, set = NULL,
                                   prefix = 'oai_dc', token = NULL, ...) {
  oai::list_identifiers(dc_oai_base(), from = from, until = until, set = set, prefix = prefix, ...)
}
