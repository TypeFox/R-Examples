#' Headers
#'
#' @export
#' @param .data Result of a call to \code{api}
#' @param ...	Key value pairs of headers
#' @param .dots	Used to work around non-standard evaluation
#' @examples \dontrun{
#' api('https://api.github.com/') %>%
#'    api_config(verbose()) %>%
#'    api_headers(`X-FARGO-SEASON` = 3) %>%
#'    peep
#'
#' api('http://httpbin.org/headers') %>%
#'    api_headers(`X-FARGO-SEASON` = 3, `X-NARCOS-SEASON` = 5)
#' }
api_headers <- function(.data, ..., .dots) {
  api_headers_(.data, .dots = lazyeval::lazy_dots(...))
}

#' @export
#' @rdname api_headers
api_headers_ <- function(.data, ..., .dots) {
  pipe_autoexec(toggle = TRUE)
  tmp <- lazyeval::all_dots(.dots, ...)
  .data <- as.req(.data)
  modifyList(.data, list(headers = getheads(tmp)))
}

getheads <- function(x) {
  tmp <- as.list(sapply(x, function(z) as.character(z$expr)))
  do.call(add_headers, tmp)
}
