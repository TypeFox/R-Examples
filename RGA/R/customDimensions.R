#' @template get_customDimensions
#' @include mgmt.R
#' @export
get_custom_dimension <- function(accountId, webPropertyId, customDimensionId, token) {
    path <- sprintf("management/accounts/%s/webproperties/%s/customDimensions/%s",
                    accountId, webPropertyId, customDimensionId)
    get_mgmt(path, token)
}

#' @template list_customDimensions
#' @include mgmt.R
#' @export
list_custom_dimensions <- function(accountId, webPropertyId, start.index = NULL, max.results = NULL, token) {
    path <- sprintf("management/accounts/%s/webproperties/%s/customDimensions",
                    accountId, webPropertyId)
    list_mgmt(path, list(start.index = start.index, max.results = max.results), token)
}
