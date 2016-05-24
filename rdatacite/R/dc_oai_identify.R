#' Learn about the Dryad OAI-PMH service
#'
#' @export
#'
#' @inheritParams dc_oai_listmetadataformats
#' @examples \dontrun{
#' dc_oai_identify()
#' }
dc_oai_identify <- function(...) {
  oai::id(dc_oai_base(), ...)
}
