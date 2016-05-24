#' Return a list of an item's parts.
#'
#' @export
#'
#' @param itemid the item id (character)
#' @inheritParams bhl_getcollections
#' @examples \dontrun{
#' bhl_getitemparts(35600)
#' bhl_getitemparts(35600, 'json')
#' bhl_getitemparts(35600, 'xml')
#' }

bhl_getitemparts <- function(itemid, as = "table", key = NULL, ...) {
  args <- bhlc(list(op = "GetItemParts", apikey = check_key(key), format = as_f(as),
                       itemid=itemid))
  bhl_GET(as, args, ...)
}
