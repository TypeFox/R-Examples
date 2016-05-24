#' Return a citation for a title, using the BibTeX format.
#'
#' @export
#' @param titleid the identifier of an individual title (numeric)
#' @inheritParams bhl_getcollections
#' @examples \dontrun{
#' bhl_gettitlebibTex(1726)
#' bhl_gettitlebibTex(1726, 'xml')
#' bhl_gettitlebibTex(1726, 'json')
#' }

bhl_gettitlebibTex <- function(titleid = NULL, as = "list", key = NULL, ...)
{
  args <- bhlc(list(op = "GetTitleBibTex", apikey = check_key(key), format = as_f(as), titleid=titleid))
  bhl_GET(as, args, ...)
}
