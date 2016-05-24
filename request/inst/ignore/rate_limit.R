#' Rate limiting
#'
#' FIXME: still need to implement doing this in the request
#'
#' @export
#'
#' @param .data Result of a call to \code{api}
#' @param value (integer) Value of rate limit, number of requests allowed
#' @param period Time period, e.g., 1 min, 60 min, 1 hr, 24 hrs
#' @param on_limit What to do on reaching rate limit. See Details.
#'
#' @details
#' \code{on_limit} options:
#' \itemize{
#'  \item stop - and give error message to use
#'  \item warn - and give error message to use
#'  \item wait - and give max time to wait
#' }
#'
#' @examples \dontrun{
#' api('https://api.github.com/') %>%
#'  api_path(repos, ropensci, rgbif, issues) %>%
#'  rate_limit(value = 5, period = "24 hrs")
#'
#' qr %>% rate_limit(value = 5, period = "24 hrs")
#' qr %>% rate_limit(value = 5000, period = "24 hrs")
#' qr %>% rate_limit(value = 10, period = "5 min")
#' qr %>% rate_limit(value = 10, period = "5 min", on_limit = with_wait(1))
#' }
rate_limit <- function(.data, value, period, on_limit = with_stop()) {
  pipe_autoexec(toggle = TRUE)
  .data <- as.req(.data)
  modifyList(.data, list(
    rate_limit = list(value = value, period = period, on_limit = on_limit))
  )
}

with_stop <- function(x = "Rate limit reached") {
  list(x = x, fxn = function(x) stop(x, call. = FALSE))
}

with_warn <- function(x = "Rate limit reached") {
  list(x = x, fxn = function(x) warning(x, call. = FALSE))
}

with_wait <- function(x = 3) {
  list(x = x, fxn = function(x) Sys.sleep(x))
}
