#' Get title metadata
#'
#' Return metadata about a title. You may choose to include a list of
#' the items (books) associated  with the title.
#'
#' @export
#' @param titleid the identifier of an individual title (numeric)
#' @param items (logical) TRUE of FALSE (default) to inclue items
#' @inheritParams bhl_getcollections
#' @examples \dontrun{
#' bhl_gettitlemetadata(1726, TRUE)
#' bhl_gettitlemetadata(1726, as='list')
#' bhl_gettitlemetadata(1726, as='xml')
#' }

bhl_gettitlemetadata <- function(titleid = NA, items = FALSE, as = "list", key = NULL, ...) {
  as <- match.arg(as, c("list","json","xml"))
  args <- bhlc(list(op = "GetTitleMetadata", apikey = check_key(key), format = as_f(as),
                       titleid = titleid, items = items))
  bhl_GET(as, args, ...)
}
