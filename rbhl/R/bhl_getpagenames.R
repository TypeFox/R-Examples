#' Return a list of names that appear on a page.
#'
#' @export
#' @param page page number to get
#' @inheritParams bhl_getcollections
#' @examples \dontrun{
#' bhl_getpagenames('1328690')
#' bhl_getpagenames('1328690', 'json')
#' bhl_getpagenames('1328690', 'list')
#' }

bhl_getpagenames <- function(page = NULL, as = 'table', key = NULL, ...)
{
  args <- bhlc(list(op = "GetPageNames", apikey=check_key(key), format=as_f(as), pageid=page))
  bhl_GET(as, args, ...)
}
