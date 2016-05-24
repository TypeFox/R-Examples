#' Private messages
#'
#' @export
#' @param username A username. required.
#' @template args
#' @examples \dontrun{
#' private_messages(username = "sckott")
#' }
private_messages <- function(username, url = NULL, key = NULL, user = NULL, ...) {
  args <- dc(list(api_key = check_key(key), api_username = check_user(user)))
  disc_GET(check_url(url), sprintf("topics/private-messages/%s.json", username), args, ...)
}
