#' @template list_accounts
#' @include mgmt.R
#' @export
list_accounts = function(start.index = NULL, max.results = NULL, token) {
    path <- "management/accounts"
    list_mgmt(path, list(start.index = start.index, max.results = max.results), token)
}
