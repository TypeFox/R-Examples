#' Notifications data
#'
#' @template args
#' @examples \dontrun{
#' notifications()
#' }
notifications <- function(url = NULL, key = NULL, user = NULL, ...) {
  args <- dc(list(api_key = check_key(key), api_username = check_user(user)))
  disc_GET(check_url(url), "notifications.json", args, ...)
}
