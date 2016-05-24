#' Use cat or dog pictures for various HTTP status codes
#'
#' @export
#' @name cat_for_status
#'
#' @param code An http status code
#' @param browse (logical) If \code{TRUE}, opens image in default browser. If
#' \code{FALSE}, returns the URL of the image.
#' @return Opens image in your default browser, or returns URL
#' @examples \dontrun{
#' # give back url
#' cat_for_status(100)
#' dog_for_status(100)
#' cat_for_status(301)
#' dog_for_status(301)
#' cat_for_status(400)
#' dog_for_status(400)
#'
#' # open image in default browser
#' cat_for_status(400, browse=TRUE)
#'
#' # not found
#' # cat_for_status(555)
#' }
cat_for_status <- function(code, browse = FALSE){
  code <- as.character(code)
  if (code %in% names(status_codes)) {
    url <- sprintf('http://httpcats.herokuapp.com/%s', code)
    if (browse) browseURL(url) else url
  } else {
    stopcode("No cat code found", "")
  }
}

#' @export
#' @rdname cat_for_status
dog_for_status <- function(code, browse = FALSE) {
  code <- as.character(code)
  if (code %in% names(status_codes)) {
    url <- sprintf('http://httpstatusdogs.com/wp-content/uploads/%s.jpg', code)
    if (browse) browseURL(url) else url
  } else {
    stopcode("No dog code found", "")
  }
}
