#' @title Retrieve the set structure of DataCite
#'
#' @description Useful for selective harvesting
#'
#' @export
#' @param token a token previously provided by the server to resume a request
#'     where it last left off.
#' @inheritParams dc_oai_listmetadataformats
#' @examples \dontrun{
#' dc_oai_listsets()
#' }
dc_oai_listsets <- function(token = NULL, ...) {
  oai::list_sets(dc_oai_base(), token = token, ...)
}
