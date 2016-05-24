#' Download metadata for individual Dryad id's
#'
#' @export
#' @param ids Dryad identifier, i.e. oai:datadryad.org:10255/dryad.8820
#' @param prefix A character string to specify the metadata format in OAI-PMH requests issued to
#' the repository. The default (\code{"oai_dc"}) corresponds to the mandatory OAI unqualified
#' Dublin Core metadata schema.
#' @param as (character) What to return. One of "df" (for data.frame; default),
#' "list", or "raw" (raw text)
#' @param ... Curl debugging options passed on to \code{\link[httr]{GET}}
#' @return XML character string, data.frame, or list, depending on what requested
#' witht the \code{as} parameter
#' @examples \dontrun{
#' dr_get_records(ids = 'oai:datadryad.org:10255/dryad.8820')
#' handles <- c('10255/dryad.36217', '10255/dryad.86943', '10255/dryad.84720', '10255/dryad.34100')
#' ids <- paste0('oai:datadryad.org:', handles)
#' dr_get_records(ids)
#' }
dr_get_records <- function(ids, prefix = "oai_dc", as = "df", ...) {
  oai::get_records(ids, prefix = prefix, url = dr_base_oai(), as = as, ...)
}
