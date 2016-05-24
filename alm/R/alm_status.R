#' Check status of an ALM service.
#' 
#' @export
#' @param key (character) Your API key, either enter, or loads from .Rprofile. Only required for 
#'    PKP source, not the others.
#' @param url API endpoint, defaults to http://alm.plos.org/api/v5/status (character)
#' @param ... optional additional curl options (debugging tools mostly) 
#' @examples \dontrun{
#' alm_status()
#' }

alm_status <- function(key = NULL, url = 'http://alm.plos.org/api/v5/status', ...)
{	
  out <- GET(url, query=almcompact(list(api_key = key)), ...)
  stop_for_status(out)
  tt <- content(out, as = "text")
  jsonlite::fromJSON(tt, TRUE)$data
}
