#' Get available Dryad metadata formats
#'
#' @export
#' @param ... Curl debugging options passed on to \code{\link[httr]{GET}}
#' @return List of information on metadata formats.
#' @examples \dontrun{
#' dr_list_metadata_formats()
#' }
dr_list_metadata_formats <- function(...) {
  oai::list_metadataformats(url = dr_base_oai(), ...)
}
