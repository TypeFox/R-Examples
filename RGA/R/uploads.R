#' @template get_uploads
#' @include mgmt.R
#' @export
get_upload <- function(accountId, webPropertyId, customDataSourceId, uploadId, token) {
    path <- sprintf("management/accounts/%s/webproperties/%s/customDataSources/%s/uploads/%s",
                    accountId, webPropertyId, customDataSourceId, uploadId)
    get_mgmt(path, token)
}

#' @template list_uploads
#' @include mgmt.R
#' @export
list_uploads <- function(accountId, webPropertyId, customDataSourceId, start.index = NULL, max.results = NULL, token) {
    path <- sprintf("management/accounts/%s/webproperties/%s/customDataSources/%s/uploads",
                    accountId, webPropertyId, customDataSourceId)
    list_mgmt(path, list(start.index = start.index, max.results = max.results), token)
}
