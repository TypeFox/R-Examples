#' Get a url, with sensible defaults
#'
#' @export
#' @param .data A request object
#' @param ... Curl options passed on to \code{\link[httr]{GET}}.
#' @param parse (logical) Attempt to parse data to data.frame if possible. Default: TRUE
#' @details Attempts to simplify the http request process by using sensible defaults:
#' \itemize{
#'  \item GET by default: you most likely want to use \code{\link[httr]{GET}}
#'  \item You most likely want a data.frame back, so we attempt to coerce to a data.frame
#' }
#' @examples \dontrun{
#' "https://api.github.com/" %>%
#'    Get()
#'
#' "https://api.github.com/" %>%
#'    Progress() %>%
#'    Verbose() %>%
#'    Get()
#'
#' "https://api.github.com/" %>%
#'    Timeout(3) %>%
#'    Get()
#'
#' "http://api.crossref.org/works/" %>%
#'    User_agent("howdydoodie") %>%
#'    Get()
#'
#' "http://api.plos.org/search" %>%
#'    Query(q = "*:*", wt = "json") %>%
#'    Get() %>%
#'    .$response %>%
#'    .$docs
#' }

Get <- function(.data, parse = TRUE, ...) {
  .data <- as.req(.data)
  .data$config <- c(httr::user_agent(make_ua()), .data$config)
  # .data$config <- c(user_agent(make_ua()), combconfig(.data$config))
  .data$url <- gather_paths(.data)
  res <- suppressWarnings(httr::GET(.data$url[1], .data$config, query = .data$query, ...))

#   hu <- httr:::handle_url(NULL, .data$url[[1]], query = .data$query)
#   req <- httr:::request_build("GET", hu$url, .data$config)
#   res <- suppressWarnings(httr:::request_perform(req, hu$handle$handle))

  # fix me, replace with error handler from .data
  if (is.null(.data$error)) {
    httr::stop_for_status(res)
  } else {
    .data$error[[1]](res)
  }
  if (grepl("json", res$headers$`content-type`)) {
    txt <- httr::content(res, "text")
    jsonlite::fromJSON(txt, parse, flatten = TRUE)
  } else {
    content(res)
  }
}

Put <- function(.data, ...) {
  .data <- as.req(.data)
  res <- httr::PUT(.data$url, body = .data$body, ...)
  httr::stop_for_status(res)
  if (grepl("json", res$headers$`content-type`)) {
    jsonlite::fromJSON(httr::content(res, "text"))
  } else {
    httr::content(res)
  }
}

#' @export
print.snapdf <- function(x, ..., n = 10){
  trunc_mat(x, n = n)
}
