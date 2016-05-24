#' Return metadata about an item.
#'
#' You may choose to include a list of the item's pages.
#'
#' @export
#'
#' @param itemid item id (character)
#' @param pages return the items pages (\code{TRUE}/\code{FALSE})
#' @param ocr (logical) \code{TRUE} to return the ocr for the item's pages. Setting this
#' to \code{TRUE} apparently doesn't return any actual ocr text, but leaving parameter
#' here for now.
#' @param parts (logical) \code{TRUE} to return the item's parts. Setting this
#' to \code{TRUE} apparently doesn't return any parts text, but leaving parameter
#' here for now.
#' @inheritParams bhl_getcollections
#'
#' @examples \dontrun{
#' bhl_getitemmetadata('16800', TRUE)
#' bhl_getitemmetadata('16800', TRUE, as='xml')
#' bhl_getitemmetadata('16800', TRUE, as='json')
#' bhl_getitemmetadata('16800', TRUE, as='list')
#'
#' # bhl_getitemmetadata(20419, pages=FALSE, parts=TRUE)
#' }

bhl_getitemmetadata <- function(itemid = NULL, pages = TRUE, ocr = FALSE, parts = FALSE,
  as = 'table', key = NULL, ...) {

  args <- bhlc(list(op = "GetItemMetadata", apikey = check_key(key), pages = pages, itemid = itemid,
                       format = as_f(as), ocr = if (ocr) 't' else NULL,
                       parts = if (parts) 't' else NULL))
  bhl_GET(as, args, ...)
}
