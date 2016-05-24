#' List identifiers of the Pangaea repository
#'
#' @export
#' @inheritParams pg_list_records
#' @param token	(character) a token previously provided by the server to resume a
#' request where it last left off. 50 is max number of records returned. We will
#' loop for you internally to get all the records you asked for.
#' @param as (character) What to return. One of "df" (for data.frame; default),
#' "list", or "raw" (raw text)
#' @return XML character string, data.frame, or list, depending on what requested
#' witht the \code{as} parameter
#' @examples \dontrun{
#' pg_list_identifiers(from='2015-09-01', until='2015-09-05')
#' pg_list_identifiers(set="geocode1", from='2015-09-01', until='2015-09-05')
#' pg_list_identifiers(prefix="iso19139", from='2015-09-01', until='2015-09-20')
#' pg_list_identifiers(prefix="dif", from='2015-09-01', until='2015-09-20')
#'
#' library('httr')
#' pg_list_identifiers(prefix="dif", from='2015-09-01', until='2015-09-05', config=verbose())
#' }

pg_list_identifiers <- function(prefix = "oai_dc", from = NULL, until = NULL,
                                set = NULL, token = NULL, as = "df", ...) {

  oai::list_identifiers(url = baseoai(), prefix = prefix, from = from, until = until,
                        set = set, token = token, as = as, ...)
}
