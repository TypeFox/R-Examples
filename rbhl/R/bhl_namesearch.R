#' Search for a particular name.
#'
#' Names both with and without NameBank identifiers are returned.
#'
#' @export
#' @param name species name (character)
#' @inheritParams bhl_getcollections
#'
#' @examples \dontrun{
#' bhl_namesearch('poa annua')
#' bhl_namesearch(name='helianthus annuus')
#' bhl_namesearch(name='helianthus annuus', as='xml')
#' bhl_namesearch(name='helianthus annuus', as='json')
#' }

bhl_namesearch <- function(name = NULL, as = "table", key = NULL, ...) {
  args <- bhlc(list(op = "NameSearch", name = name, apikey = check_key(key), format = as_f(as)))
  bhl_GET(as, args, ...)
}
