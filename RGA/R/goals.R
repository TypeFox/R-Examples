#' @template get_goals
#' @include mgmt.R
#' @export
get_goal <- function(accountId, webPropertyId, profileId, goalId, token) {
    path <- sprintf("management/accounts/%s/webproperties/%s/profiles/%s/goals/%s",
                    accountId, webPropertyId, profileId, goalId)
    get_mgmt(path, token)
}

#' @template list_goals
#' @include mgmt.R
#' @export
list_goals = function(accountId = "~all", webPropertyId = "~all", profileId = "~all", start.index = NULL, max.results = NULL, token) {
    path <- sprintf("management/accounts/%s/webproperties/%s/profiles/%s/goals",
                    accountId, webPropertyId, profileId)
    list_mgmt(path, list(start.index = start.index, max.results = max.results), token)
}
