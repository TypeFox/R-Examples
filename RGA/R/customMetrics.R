#' @template get_customMetrics
#' @include mgmt.R
#' @export
get_custom_metric <- function(accountId, webPropertyId, customMetricId, token) {
    path <- sprintf("management/accounts/%s/webproperties/%s/customMetrics/%s",
                    accountId, webPropertyId, customMetricId)
    get_mgmt(path, token)
}

#' @template list_customMetrics
#' @include mgmt.R
#' @export
list_custom_metrics <- function(accountId, webPropertyId, start.index = NULL, max.results = NULL, token) {
    path <- sprintf("management/accounts/%s/webproperties/%s/customMetrics",
                    accountId, webPropertyId)
    list_mgmt(path, list(start.index = start.index, max.results = max.results), token)
}
