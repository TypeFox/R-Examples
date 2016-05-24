#' Listing BigML Sources
#' @export
#' @family source methods
#' @references \url{https://bigml.com/developers/sources}
#' @param flatten A logical value indicating whether to flatten the response
#'	into a dataframe.
#' @param sources_only A logical value indicating whether to only return
#'	the data frame of source information (only valid if \code{flatten} is
#'	\code{TRUE}).
#' @template dots
#' @template author
#' @template source_list_return
listSources <-
function (flatten = TRUE, sources_only = TRUE, ...)
{
    message("Retrieving the source list...")
    response = .basic_api(.SOURCE_URL)$list(...)
    if (flatten) {
        response$sources = ldply(response$objects, function(x) {
            x$fields = NULL
            x$source_parser$missing_tokens = paste(x$source_parser$missing_tokens,
                collapse = ",")
            as.data.frame(unlist(x, recursive = F), stringsAsFactors=FALSE)
        })
        response$fields = ldply(response$objects, function(y) {
            ldply(y$fields, function(z) {
                z$resource = y$resource
                data.frame(z, stringsAsFactors=FALSE)
            })
        })
        response$objects = NULL
        if (sources_only) {
            response = response$sources
        }
    }
    return(response)
}
