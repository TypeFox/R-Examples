#' Get a list of languages in which materials in BHL have been written.
#'
#' @export
#' @inheritParams bhl_getcollections
#' @examples \dontrun{
#' bhl_getlanguages()
#' bhl_getlanguages('json')
#' bhl_getlanguages('xml')
#' bhl_getlanguages('list')
#' }
bhl_getlanguages <- function(as = "table", key = NULL, ...)
{
  args <- bhlc(list(op = "GetLanguages", apikey = check_key(key), format = as_f(as)))
  bhl_GET(as, args, ...)
}
