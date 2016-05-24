#' Return a list of the identifiers of all unpublished items.
#'
#' @export
#'
#' @param type The type of identifier (doi, oclc, issn, isbn, lccn, ddc, nal, nlm, coden)
#' @param value The identifier value
#' @inheritParams bhl_getcollections
#' @examples \dontrun{
#' bhl_getpartbyidentifier('doi', '10.4039/Ent38406-12')
#' bhl_getpartbyidentifier('doi', '10.4039/Ent38406-12', as='json')
#' bhl_getpartbyidentifier('doi', '10.4039/Ent38406-12', as='xml')
#' bhl_getpartbyidentifier('doi', '10.4039/Ent38406-12', as='list')
#' }

bhl_getpartbyidentifier <- function(type=NULL, value=NULL, as = "table", key = NULL, ...)
{
  args <- bhlc(list(op = "GetPartByIdentifier", apikey = check_key(key), format = as_f(as),
                       type=type, value=value))
  bhl_GET(as, args, ...)
}
