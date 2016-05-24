#' @template list_segments
#' @include mgmt.R
#' @export
list_segments = function(start.index = NULL, max.results = NULL, token) {
    path <- "management/segments"
    list_mgmt(path, list(start.index = start.index, max.results = max.results), token)
}
