#' Learn about the Dryad OAI-PMH service.
#'
#' @export
#' @param ... Curl debugging options passed on to \code{\link[httr]{GET}}
#' @return List of information describing Dryad.
#' @examples \dontrun{
#' dr_identify()
#' }
dr_identify <- function(...) {
  oai::id(url = dr_base_oai(), ...)
}
