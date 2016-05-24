#' Find and return metadata about an item or items that match a specific identifier.
#'
#' If you know the Internet Archive identifier for an item, use this method to
#'    look up the equivalent item in BHL.
#'
#' @export
#'
#' @param type the type of identifier (barcode or ia) (character)
#' @param value the identifier value (character)
#' @inheritParams bhl_getcollections
#'
#' @examples \dontrun{
#' bhl_getitembyidentifier(type='ia', value='animalkingdomarr03cuvi')
#' bhl_getitembyidentifier(type='ia', value='animalkingdomarr03cuvi', as='json')
#' bhl_getitembyidentifier(type='ia', value='animalkingdomarr03cuvi', as='xml')
#' }

bhl_getitembyidentifier <- function(type = NULL, value = NULL, as = 'table',
  key = NULL, ...) {

  args <- bhlc(list(op = "GetItemByIdentifier", apikey = check_key(key), type = type,
                       value = value, format = as_f(as)))
  bhl_GET(as, args, ...)
}
