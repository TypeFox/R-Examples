#' http DSL
#'
#' @export
#' @family dsl
#' @examples \dontrun{
#' # Set base url
#' ## works with full or partial URLs
#' api('https://api.github.com/')
#' api('http://api.gbif.org/v1')
#' api('api.gbif.org/v1')
#'
#' ## works with ports, full or partial
#' api('http://localhost:9200')
#' api('localhost:9200')
#' api(':9200')
#' api('9200')
#' api('9200/stuff')
#'
#' # set paths
#' ## NSE
#' api('https://api.github.com/') %>%
#'   api_path(repos, ropensci, rgbif, issues)
#' ## SE
#' api('https://api.github.com/') %>%
#'   api_path_('repos', 'ropensci', 'rgbif', 'issues')
#'
#' # template
#' repo_info <- list(username = 'craigcitro', repo = 'r-travis')
#' api('https://api.github.com/') %>%
#'   api_template(template = 'repos/{{username}}/{{repo}}/issues', data = repo_info)
#'
#' # OAuth setup
#' api('https://api.github.com/') %>%
#'  api_oauth2(token = "<token>")
#'
#' # Error handler
#' api('https://api.github.com/') %>%
#'  api_error_handler(stop_for_status)
#'
#' # Config handler
#' api('https://api.github.com/') %>%
#'  api_config(verbose())
#'
#' # Query handler
#' ## NSE
#' api("http://api.plos.org/search") %>%
#'   api_query(q = ecology, wt = json, fl = 'id,journal') %>%
#'   Get()
#' ## SE
#' api("http://api.plos.org/search") %>%
#'   api_query_(q = "ecology", wt = "json", fl = 'id', fl = 'journal') %>%
#'   Get()
#'
#' # xxx handler
#' api("http://api.plos.org/search") %>%
#'   xxxx
#'
#' # Full examples
#' api('https://api.github.com/') %>%
#'   api_path(repos, ropensci, rgbif, issues) %>%
#'   api_config(verbose()) %>%
#'   Get()
#'
#' repo_info <- list(username = 'craigcitro', repo = 'r-travis')
#' api('https://api.github.com/') %>%
#'   api_template(template = 'repos/{{username}}/{{repo}}/issues', data = repo_info) %>%
#'   Get()
#'
#' # parse=TRUE by default, in this eg parses directly to a data.frame
#' api('https://api.github.com/') %>%
#'   api_path(repos, ropensci, rgbif, issues) %>%
#'   Get()
#' }
api <- function(x) {
  structure(list(url = as.url(x)), class = "endpoint")
}

print.endpoint <- function(x, ...) {
  cat(sprintf("URL: %s", x$url))
}

#' path defintion ------------------------------------
#' @export
#' @family dsl
api_path <- function(.data, ..., .dots) {
  api_path_(.data, .dots = lazyeval::lazy_dots(...))
}

#' @export
#' @family dsl
api_path_ <- function(.data, ..., .dots) {
  tmp <- lazyeval::all_dots(.dots, ...)
  .data <- as.req(.data)
  modifyList(.data, list(paths = getpaths(tmp)))
}

getpaths <- function(x) {
  unname(sapply(x, function(z) as.character(z$expr)))
}

#' api template ------------------------------------
#' @export
#' @family dsl
api_template <- function(.data, template, data) {
  .data <- as.req(.data)
  temp <- whisker::whisker.render(template, data)
  modifyList(.data, list(template = temp))
}

#' oauth setup ------------------------------------
#' @export
#' @family dsl
api_oauth2 <- function(.data, token = NULL, app_name = NULL, key = NULL,
                       secret = NULL, request = NULL, authorize = NULL,
                       access = NULL, base_url = NULL, ...) {
  .data <- as.req(.data)

  args <- comp(list(token = token, app_name = app_name, key = key, secret = secret))
  if (length(args) == 0) {
    stop("either token or app_name + key + secret must be provided", call. = FALSE)
  } else {
    if (!is.null(token)) {
      auth <- config(token = token)
    } else {
      app <- oauth_app(app_name, key, secret)
      endpts <- oauth_endpoint(request = request, authorize = authorize,
                               access = access, base_url = base_url)
      token <- oauth2.0_token(endpts, app)
      auth <- config(token = token)
    }
  }

  modifyList(.data, list(config = c(auth)))
}

#' error handler ------------------------------------
#' @export
#' @family dsl
api_error_handler <- function(.data, func) {
  .data <- as.req(.data)
  fn_name <- deparse(substitute(func))
  tmp <- setNames(list(func), fn_name)
  modifyList(.data, list(error = tmp))
}

#' configuration handler ------------------------------------
#' @export
#' @family dsl
api_config <- function(.data, ...) {
  .data <- as.req(.data)
  tmp <- list(...)
  modifyList(.data, list(config = tmp))
}

#' query handler ------------------------------------
#' @export
#' @family dsl
api_query <- function(.data, ...){
  api_query_(.data, .dots = lazyeval::lazy_dots(...))
}

api_query_ <- function(.data, ..., .dots){
  dots <- lazyeval::all_dots(.dots, ...)
  args <- sapply(dots, "[[", "expr")
  .data <- as.req(.data)
  modifyList(.data, list(query = args))
}

# api_query <- function(.data, ...){
#   .data <- as.req(.data)
#   args <- list(...)
#   modifyList(.data, list(query = args))
# }

# #' @export
# api_config <- function(.data, ..., .dots) {
#   api_config_(.data, .dots = lazyeval::lazy_dots(...))
# }
#
# #' @export
# api_config_ <- function(.data, ..., .dots) {
#   tmp <- lazyeval::all_dots(.dots, ...)
#   .data <- as.req(.data)
#   modifyList(.data, list(config = tmp))
# }
