#' Return a list of the identifiers of all unpublished items.
#'
#' @export
#' @inheritParams bhl_getcollections
#' @examples \dontrun{
#' bhl_getunpublisheditems()
#' bhl_getunpublisheditems('xml')
#' bhl_getunpublisheditems('json')
#' }

bhl_getunpublisheditems <- function(as = "table", key = NULL, ...)
{
  args <- bhlc(list(op = "GetUnpublishedItems", apikey = check_key(key), format = as_f(as)))
  bhl_GET(as, args, ...)
}
