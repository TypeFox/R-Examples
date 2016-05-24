#' Email
#'
#' @name email
#' @export
#' @param filter An email filter, one of received, sent, skipped, rejected
#' @template args
#' @examples \dontrun{
#' # email settings
#' email_settings()
#'
#' # list settings for an email address
#' list_email("received")
#' }
email_settings <- function(url = NULL, key = NULL, user = NULL, ...) {
  args <- dc(list(api_key = check_key(key), api_username = check_user(user)))
  disc_GET(check_url(url), "admin/email.json", args, ...)
}

#' @export
#' @rdname email
list_email <- function(filter, url = NULL, key = NULL, user = NULL, ...) {
  args <- dc(list(api_key = check_key(key), api_username = check_user(user)))
  disc_GET(check_url(url), sprintf("admin/email/%s.json", filter), args, ...)
}
