#' List available metadata formats that the DataCite repository can disseminate
#'
#' @export
#' @param id DataCite identifier, e.g., "oai:oai.datacite.org:6718729". If left
#' blank, get all metadataformats
#' @param ... Curl options passed on to \code{\link[httr]{GET}}
#' @examples \dontrun{
#' dc_oai_listmetadataformats()
#' dc_oai_listmetadataformats("oai:oai.datacite.org:6718729")
#' }
dc_oai_listmetadataformats <- function(id = NULL, ...) {
  oai::list_metadataformats(dc_oai_base(), id, ...)
}
