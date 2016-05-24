#' Return a citation for a part, using the EndNote format.
#'
#' @export
#' @param partid The identifier of an individual part (article, chapter, etc) (numeric)
#' @inheritParams bhl_getcollections
#' @examples \dontrun{
#' bhl_getpartendnote(1000)
#' bhl_getpartendnote(1000, as='xml')
#' bhl_getpartendnote(1000, as='json')
#' }

bhl_getpartendnote <- function(partid, as = "list", key = NULL, ...)
{
  args <- bhlc(list(op = "GetPartEndNote", apikey = check_key(key), format = as_f(as), partid=partid))
  bhl_GET(as, args, ...)
}
