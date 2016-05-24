#' Retrieving a BigML Source
#' @export
#' @family source methods
#' @references \url{https://bigml.com/developers/sources}
#' @param source_id A character value giving the name of the source.
#' @param flatten A logical value indicating whether to flatten the response
#'	into a data frame.
#' @template dots
#' @template author
#' @template normal_methods
#' @template source_return
getSource <-
function (source_id, flatten = TRUE)
{
    message("Retrieving the source...")
    response = .basic_api(.SOURCE_URL)$get(source_id)
    if (flatten) {
        response$fields = ldply(response$fields, function(x) {
            data.frame(x, stringsAsFactors=FALSE)
        })
    }
    return(response)
}
