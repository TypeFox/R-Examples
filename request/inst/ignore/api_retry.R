#' Retry on failed request N times
#'
#' FIXME: still need to implement doing this in the request
#'
#' @export
#' @param .data Result of a call to \code{api}
#' @param n Number of times to repeat the request
#' @param time Number of seconds to delay between repeated calls
#'
#' @details This doesn't use the retry option within curl itself,
#' as it's not available via the \code{curl} R client. Instead, we
#' retry X times you specify, if the previous call failed.
#'
#' @examples \dontrun{
#' api('http://127.0.0.1:8080') %>%
#'    api_path(repos, asdfasdf) %>%
#'    api_retry(n = 5) %>%
#'    http
#'
#' api('http://127.0.0.1:8080') %>%
#'    api_retry(n = 5, time = 2) %>%
#'    http
#' }
api_retry <- function(.data, n, time = 1) {
  .data <- as.req(.data)
  modifyList(.data, list(retry = list(n = n, time = time)))
}
