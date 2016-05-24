#' @template get_webproperties
#' @include mgmt.R
#' @export
get_webproperty <- function(accountId, webPropertyId, token) {
    path <- sprintf("management/accounts/%s/webproperties/%s",
                    accountId, webPropertyId)
    get_mgmt(path, token)
}

#' @template list_webproperties
#' @include mgmt.R
#' @export
list_webproperties = function(accountId = "~all", start.index = NULL, max.results = NULL, token) {
    path <- sprintf("management/accounts/%s/webproperties", accountId)
    list_mgmt(path, list(start.index = start.index, max.results = max.results), token)
}
