#' Tag data
#'
#' @param tag A tag name
#' @template args
#' @examples \dontrun{
#' tag('r')
#' tag('geospatial')
#' }
tag <- function(tag, url = NULL, key = NULL, user = NULL, ...) {
  args <- dc(list(api_key = check_key(key), api_username = check_user(user)))
  disc_GET(check_url(url), sprintf("tags/%s.json", tag), args, ...)
}
