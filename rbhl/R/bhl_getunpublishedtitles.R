#' Return a list of the identifiers of all unpublished titles.
#'
#' @export
#' @inheritParams bhl_getcollections
#' @examples \dontrun{
#' bhl_getunpublishedtitles()
#' bhl_getunpublishedtitles('json')
#' bhl_getunpublishedtitles('xml')
#' }

bhl_getunpublishedtitles <- function(as = "table", key = NULL, ...)
{
  args <- bhlc(list(op = "GetUnpublishedTitles", apikey = check_key(key), format = as_f(as)))
  bhl_GET(as, args, ...)
}
