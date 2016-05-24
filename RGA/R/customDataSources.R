#' @template list_customDataSources
#' @include mgmt.R
#' @export
list_custom_data_sources <- function(accountId, webPropertyId, start.index = NULL, max.results = NULL, token) {
    path <- sprintf("management/accounts/%s/webproperties/%s/customDataSources",
                    accountId, webPropertyId)
    list_mgmt(path, list(start.index = start.index, max.results = max.results), token)
}
