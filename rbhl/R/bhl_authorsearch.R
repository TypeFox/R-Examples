#' Return a list of authors that match (fully or partially) the specified
#'    search string.
#'
#' The namessearched are those contained in MARC 100a, 110a, 111a, 700a,
#'    710a, and 711a library records.
#'
#' @export
#' @inheritParams bhl_getcollections
#' @param name full or partial name of the author for which to search
#'     (last name listed first, i.e. 'Darwin, Charles') (character)
#' @examples \dontrun{
#' bhl_authorsearch(name='dimmock')
#' bhl_authorsearch(name='Jones')
#' }

bhl_authorsearch <- function(name = NULL, as='table', key = NULL, ...) {
  args <- bhlc(list(op = "AuthorSearch", name = name, apikey = check_key(key), format = as_f(as)))
  bhl_GET(as, args, ...)
}
