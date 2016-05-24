#' Get record from the Pangaea repository
#'
#' @export
#' @param identifier Dataset identifier. See Examples.
#' @param prefix A character string to specify the metadata format in OAI-PMH requests issued to
#' the repository. The default (\code{"oai_dc"}) corresponds to the mandatory OAI unqualified
#' Dublin Core metadata schema.
#' @param as (character) What to return. One of "df" (for data.frame; default),
#' "list", or "raw" (raw text)
#' @param ... Curl debugging options passed on to \code{\link[httr]{GET}}
#' @return XML character string, data.frame, or list, depending on what requested
#' witht the \code{as} parameter
#' @examples \dontrun{
#' pg_get_record(identifier = "oai:pangaea.de:doi:10.1594/PANGAEA.788382")
#' pg_get_record(identifier = "oai:pangaea.de:doi:10.1594/PANGAEA.269656",
#' prefix="iso19139")
#' pg_get_record(identifier = "oai:pangaea.de:doi:10.1594/PANGAEA.269656",
#' prefix="dif")
#'
#' # curl options
#' library('httr')
#' pg_get_record(identifier = "oai:pangaea.de:doi:10.1594/PANGAEA.788382", config=verbose())
#'
#' # invalid record id
#' # pg_get_record(identifier = "oai:pangaea.de:doi:10.1594/PANGAEA.11111")
#' # pg_get_record(identifier = "oai:pangaea.de:doi:10.1594/PANGAEA.11111", prefix="adfadf")
#' }
pg_get_record <- function(identifier, prefix = "oai_dc", as = "df", ...){
  oai::get_records(ids = identifier, prefix = prefix, url = baseoai(), as = as, ...)
}
