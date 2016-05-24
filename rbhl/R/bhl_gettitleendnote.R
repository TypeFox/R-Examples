#' Return a citation for a title, using the EndNote format.
#'
#' @export
#' @param titleid the identifier of an individual title (numeric)
#' @param key your BHL API key, either enter, or loads from .Rprofile
#' @param ... Curl options passed on to \code{\link[httr]{GET}}
#' @inheritParams bhl_getcollections
#' @examples \dontrun{
#' bhl_gettitleendNote(1726)
#' }

bhl_gettitleendNote <- function(titleid = NA, key = NULL, ...) {
  args <- bhlc(list(op = "GetTitleEndNote", apikey = check_key(key),
                    format = 'json', titleid = titleid))
  if (length(args) == 0) args <- NULL
  out <- GET(bhl_url(), query = args, ...)
  stop_for_status(out)
  tt <- content(out)
  gsub("\n|%.{1}", "", tt$Result)
}
