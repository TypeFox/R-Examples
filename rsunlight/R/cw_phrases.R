#' Capitol words phrases.json method. List the top phrases for a facet.
#'
#' @export
#' @param entity_type The entity type to get top phrases for. One of 'date',
#'    'month', 'state', or 'legislator'. Required.
#' @param entity_value The value of the entity given in \code{entity_type}.
#'    See Details. Required.
#' @param n The size of phrase, in words, to search for (up to 5).
#' @param page The page of results to show. 100 results are shown at a time. To get more than
#' 100 results, use the page parameter.
#' @param per_page Number of records to return. Default: 20. Max: 50.
#' @param sort The value on which to sort the results. You have to specify ascending or
#' descending (see details), but if you forget, we make it ascending by default
#' (prevents 500 error :)). Valid values are 'tfidf' (default) and 'count'.
#' @param key Your SunlightLabs API key; loads from .Rprofile.
#' @param ... Further curl options (debugging tools mostly)
#' @param as (character) One of table (default), list, or response (httr response object).
#' @return Data frame of observations by date.
#' @examples \dontrun{
#' cw_phrases(entity_type='month', entity_value=201007)
#' cw_phrases(entity_type='state', entity_value='NV')
#' cw_phrases(entity_type='legislator', entity_value='L000551')
#'
#' library('httr')
#' head(cw_phrases(entity_type='month', entity_value=201007, config=verbose()))
#' }
cw_phrases <- function(entity_type, entity_value, n = NULL, page = NULL,
  per_page = NULL, sort = NULL, as = 'table', key = NULL, ...) {

  key <- check_key(key)
  if (!is.null(sort)) {
    if (!grepl('asc|desc', sort)) {
      sort <- paste(sort, 'asc')
    }
  }
  args <- sc(list(apikey = key, entity_type = entity_type, entity_value = entity_value,
                          n = n, page = page, per_page = per_page, sort = sort))
  give_noiter(as, cwurl(), "/phrases.json", args, ...)
}
