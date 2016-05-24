#' Get a list of collections which are used to group titles and items. A single
#' collection may contain either titles or items, but not both.
#'
#' @export
#' @param as (character) Return a list ("list"), json ("json"), xml ("xml"), or parsed table
#' ("table", default). Note that \code{as="table"} can give different data format back
#' depending on the function - for example, sometimes a data.frame and sometimes a
#' character vector.
#' @param key Your BHL API key, either enter, or loads from your \code{.Renviron} as \code{BHL_KEY}
#' or from \code{.Rprofile} as \code{bhl_key}.
#' @param ... Curl options passed on to \code{\link[httr]{GET}}
#'
#' @examples \dontrun{
#' bhl_getcollections()
#' bhl_getcollections(as = 'list')
#' bhl_getcollections(as = 'json')
#' bhl_getcollections(as = 'xml')
#' }

bhl_getcollections <- function(as = 'table', key = NULL, ...) {
  args <- bhlc(list(op = "GetCollections", apikey = check_key(key), format = as_f(as)))
  bhl_GET(as, args, ...)
}
