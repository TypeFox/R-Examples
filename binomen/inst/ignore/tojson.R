#' Convert to JSON
#'
#' NOT WORKING YET
#'
#' @keywords internal
#' @param input Input object
#' @param ... Further args passed on to \code{\link[jsonlite]{toJSON}}
#' @examples \dontrun{
#' out <- make_taxon(genus="Poa", epithet="annua", authority="L.")
#' x <- taxonref(rank="species", name="Homo sapiens", id=3454, uri="http://things.com")
#' as.jsonld(x)
#' }

as.jsonld <- function(input, ...){
  jsonlite::toJSON(input, auto_unbox = TRUE, ...)
}
