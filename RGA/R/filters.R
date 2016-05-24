#' @template get_filters
#' @include mgmt.R
#' @export
get_filter <- function(accountId, filterId, token) {
    path <- sprintf("management/accounts/%s/filters/%s", accountId, filterId)
    get_mgmt(path, token)
}

#' @template list_filters
#' @include mgmt.R
#' @export
list_filters <- function(accountId, start.index = NULL, max.results = NULL, token) {
    path <- sprintf("management/accounts/%s/filters", accountId)
    list_mgmt(path, list(start.index = start.index, max.results = max.results), token)
}
