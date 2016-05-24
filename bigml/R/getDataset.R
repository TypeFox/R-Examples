#' Retrieving a BigML Dataset
#' @export
#' @family dataset methods
#' @references \url{https://bigml.com/developers/datasets}
#' @param source_id A string giving the name of the source id.
#' @param include_overview A logical value indicating whether to provide a
#'	simple data frame overview of fields.
#' @template dots
#' @template author
#' @template dataset_return
#' @template normal_methods
getDataset <-
function (source_id, include_overview = TRUE, ...)
{
    message("Retrieving the dataset...")
    response = .basic_api(.DATASET_URL)$get(id = source_id)
    if (include_overview) {
        tmpfields = ldply(response$fields, function(y) {
            y$summary = NULL
            as.data.frame(y, stringsAsFactors=FALSE)
        })
        response$fields_overview = tmpfields
    }
    return(response)
}
