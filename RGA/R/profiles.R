#' @template get_profiles
#' @include mgmt.R
#' @export
get_profile <- function(accountId, webPropertyId, profileId, token) {
    path <- sprintf("management/accounts/%s/webproperties/%s/profiles/%s",
                    accountId, webPropertyId, profileId)
    get_mgmt(path, token)
}

#' @template list_profiles
#' @include mgmt.R
#' @export
list_profiles = function(accountId = "~all", webPropertyId = "~all", start.index = NULL, max.results = NULL, token) {
    path <- sprintf("management/accounts/%s/webproperties/%s/profiles",
                    accountId, webPropertyId)
    list_mgmt(path, list(start.index = start.index, max.results = max.results), token)
}
