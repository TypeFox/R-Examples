#' @template get_experiments
#' @include mgmt.R
#' @export
get_experiment <- function(accountId, webPropertyId, profileId, experimentId, token) {
    path <- sprintf("management/accounts/%s/webproperties/%s/profiles/%s/experiments/%s",
                    accountId, webPropertyId, profileId, experimentId)
    get_mgmt(path, token)
}

#' @template list_experiments
#' @include mgmt.R
#' @export
list_experiments <- function(accountId, webPropertyId, profileId, start.index = NULL, max.results = NULL, token) {
    path <- sprintf("management/accounts/%s/webproperties/%s/profiles/%s/experiments",
                    accountId, webPropertyId, profileId)
    list_mgmt(path, list(start.index = start.index, max.results = max.results), token)
}
