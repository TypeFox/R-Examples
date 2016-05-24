#' Return a citation for a part, using the BibTeX format.
#'
#' @export
#' @param partid The identifier of an individual part (article, chapter, etc) (numeric)
#' @inheritParams bhl_getcollections
#' @examples \dontrun{
#' bhl_getpartbibtex(1000)
#' bhl_getpartbibtex(1000, 'xml')
#' bhl_getpartbibtex(1000, 'json')
#' }

bhl_getpartbibtex <- function(partid, as = "list", key = NULL, ...)
{
  args <- bhlc(list(op = "GetPartBibTex", apikey = check_key(key), format = as_f(as), partid=partid))
  bhl_GET(as, args, ...)
}
