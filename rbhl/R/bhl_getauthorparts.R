#' Return a list of parts (articles, chapters, etc) associated with a given BHL
#' author identifier. Unless the identifier  for a particular BHL author record
#' is known in advance, this method should be used in combination	with the
#' AuthorSearch method.
#'
#' Note: haven't seen examples for this function that work yet...
#'
#' @export
#' @param creatorid BHL identifier for a particular author (numeric)
#' @inheritParams bhl_getcollections
#'
#' @examples \dontrun{
#' # bhl_getauthorparts(147)
#' # bhl_getauthorparts(39120, as='json')
#' # bhl_getauthorparts(39120, as='xml')
#' # bhl_getauthorparts(39120, as='list')
#' }
bhl_getauthorparts <- function(creatorid, as='table', key = NULL, ...) {

  args <- bhlc(list(op = "GetAuthorParts", apikey = check_key(key), format = as_f(as),
                       creatorid = creatorid))
  bhl_GET(as, args, ...)
}
