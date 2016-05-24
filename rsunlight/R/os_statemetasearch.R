#' Search OpenStates metadata.
#'
#' @param state One or more two-letter state abbreviations (character)
#' @param key your SunlightLabs API key; or loads from .Rprofile
#' @param ... Curl options passed on to \code{\link[httr]{GET}}
#' @return A list with metadata for each state.
#' @export
#' @examples \dontrun{
#' os_statemetasearch()
#' os_statemetasearch('ca')
#' os_statemetasearch('tx')
#' out <- os_statemetasearch(c('tx','nv'))
#' out$tx
#' out$nv
#'
#' library('httr')
#' os_statemetasearch('tx', config=verbose())
#' }
os_statemetasearch <- function(state = NULL, key = NULL, ...) {
  key <- check_key(key)
  getdata <- function(x = NULL){
    url <- paste0(osurl(), "/metadata/")
    if (!is.null(x)) {
      url <- paste(url, x, sep = "")
    }
    url2 <- paste(url, "?apikey=", key, sep = "")
    tt <- GET(url2, ...)
    stop_for_status(tt)
    out <- content(tt, as = "text")
    jsonlite::fromJSON(out)
  }

  if (!is.null(state)) {
    res <- lapply(state, getdata)
    names(res) <- state
    res
  } else {
    getdata()
  }
}
