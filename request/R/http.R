#' Make a HTTP request
#'
#' @export
#'
#' @param req A \code{req} class object
#' @param method (character) Pick which HTTP method to use. Only GET and
#' POST for now. Default: GET
#'
#' @details By default, a GET request is made. Will fix this soon to easily allow
#' a different HTTP verb.
#'
#' The \code{http} function makes the request and gives back the parsed result.
#' Whereas, the \code{http_client} function makes the request, but gives back
#' the raw R6 class object, which you can inspect all parts of, modify, etc.
#' @examples \dontrun{
#' # high level - http()
#' api('https://api.github.com/') %>%
#'   api_path(repos, ropensci, rgbif, commits) %>%
#'   http()
#'
#' # low level - http_client()
#' res <- api('https://api.github.com/') %>%
#'   api_path(repos, ropensci, rgbif, commits) %>%
#'   http_client()
#' res$count()
#' res$body()
#' res$status()
#' res$result
#' res$links
#' res$parse()
#'
#' # Specify HTTP verb
#' api("http://httpbin.org/post") %>%
#'    api_body(x = "A simple text string") %>%
#'    http("POST")
#' }
http <- function(req, method = "GET") {
  pipe_autoexec(toggle = FALSE)
  if (!method %in% c("GET", "POST")) stop("method must be one of GET or POST", call. = FALSE)
  rr <- RequestIterator$new(limit = req$paging$limit, limit_max = req$paging$limit_max)
  switch(method,
         GET = rr$GET(req),
         POST = rr$POST(req)
  )
  rr$parse()
}

#' @export
#' @rdname http
http_client <- function(req) {
  pipe_autoexec(toggle = FALSE)
  rr <- RequestIterator$new(limit = req$paging$limit, limit_max = req$paging$limit_max)
  rr$GET(req)
  return(rr)
}
