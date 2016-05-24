#' Return a list of the identifiers of all unpublished parts (articles, chapters, etc).
#'
#' @export
#' @inheritParams bhl_getcollections
#' @examples \dontrun{
#' bhl_getunpublishedparts()
#' bhl_getunpublishedparts('json')
#' bhl_getunpublishedparts('xml')
#' }

bhl_getunpublishedparts <- function(as = "table", key = NULL, ...)
{
  args <- bhlc(list(op = "GetUnpublishedParts", apikey = check_key(key), format = as_f(as)))
  bhl_GET(as, args, ...)
}
