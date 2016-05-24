#' List the sets in the Dryad metadata repository.
#'
#' Retrieve the set structure of Dryad, useful for selective harvesting
#'
#' @export
#' @param token	(character) a token previously provided by the server to resume a
#' request where it last left off. 50 is max number of records returned. We will
#' loop for you internally to get all the records you asked for.
#' @param as (character) What to return. One of "df" (for data.frame; default),
#' "list", or "raw" (raw text)
#' @param ... Curl debugging options passed on to \code{\link[httr]{GET}}
#' @examples \dontrun{
#' dr_list_sets()
#' dr_list_sets(as = "list")
#' dr_list_sets(as = "raw")
#' library('httr')
#' res <- dr_list_sets(config = verbose())
#' }
dr_list_sets <- function(token = NULL, as = "df", ...) {
  oai::list_sets(url = dr_base_oai(), token = token, as = as, ...)
}
