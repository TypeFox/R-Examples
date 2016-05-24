#' Perform a simple title search.
#'
#' The full title (as specified in MARC 245a and MARC 245b library records)
#'    is searched for the specified string. Basic metadata for all full and
#'    partial matches is returned.
#'
#' @export
#'
#' @param title full or partial title for which to search (character)
#' @inheritParams bhl_getcollections
#' @examples \dontrun{
#' bhl_titlesearchsimple('nematocerous')
#' bhl_titlesearchsimple('husbandry')
#' }
bhl_titlesearchsimple <- function(title = NA, as = "table", key = NULL, ...)
{
  args <- bhlc(list(op = "TitleSearchSimple", apikey = check_key(key), format = as_f(as), title=title))
  bhl_GET(as, args, ...)
}
