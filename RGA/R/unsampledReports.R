#' @template get_unsampledReports
#' @include mgmt.R
#' @export
get_unsampled_report <- function(accountId, webPropertyId, profileId, unsampledReportId, token) {
    path <- sprintf("management/accounts/%s/webproperties/%s/profiles/%s/unsampledReports/%s",
                    accountId, webPropertyId, profileId, unsampledReportId)
    get_mgmt(path, token)
}

#' @template list_unsampledReports
#' @include mgmt.R
#' @export
list_unsampled_reports <- function(accountId, webPropertyId, profileId, start.index = NULL, max.results = NULL, token) {
    path <- sprintf("management/accounts/%s/webproperties/%s/profiles/%s/unsampledReports",
                    accountId, webPropertyId, profileId)
    list_mgmt(path, list(start.index = start.index, max.results = max.results), token)
}
