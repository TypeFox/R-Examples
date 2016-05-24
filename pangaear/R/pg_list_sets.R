#' List the set structure of the Pangaea repository
#'
#' @export
#' @param token	(character) a token previously provided by the server to resume a
#' request where it last left off. 50 is max number of records returned. We will
#' loop for you internally to get all the records you asked for.
#' @param as (character) What to return. One of "df" (for data.frame; default),
#' "list", or "raw" (raw text)
#' @param ... Curl debugging options passed on to \code{\link[httr]{GET}}
#' @return XML character string, data.frame, or list, depending on what requested
#' witht the \code{as} parameter
#' @examples \dontrun{
#' pg_list_sets()
#' pg_list_sets(as = "list")
#' pg_list_sets(as = "raw")
#' library('httr')
#' res <- pg_list_sets(config = verbose())
#' }

pg_list_sets <- function(token = NULL, as = "df", ...) {
  oai::list_sets(url = baseoai(), token = token, as = as, ...)
}
