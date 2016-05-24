#' Listing BigML Datasets
#' @export
#' @family dataset methods
#' @references \url{https://bigml.com/developers/datasets}
#' @param flatten A logical value indicating whether to flatten the response
#'	into a dataframe.
#' @param datasets_only A logical value indicating whether to only return
#'	the data frame of field information (only valid if \code{flatten} is
#'	\code{TRUE}).
#' @template dots
#' @template dataset_list_return
#' @template author
listDatasets <-
function (flatten = TRUE, datasets_only = TRUE, ...)
{
    message("Retrieving the dataset list...")
    response = .basic_api(.DATASET_URL)$list(...)
    if (flatten) {
        response$datasets = ldply(response$objects, function(x) {
            x$fields = NULL
            x$field_types = NULL
            x$field_new_fields = NULL
            x$status$field_errors = NULL
            x$status$row_format_errors = NULL
            x$output_fields = NULL
            x$tags = NULL
            print (unlist(x, recursive = F))
            as.data.frame(unlist(x, recursive = F), stringsAsFactors=FALSE)
        })
        response$fields = ldply(response$objects, function(y) {
            ldply(y$fields, function(z) {
                z = as.list(z)
                z$resource = y$resource
                as.data.frame(z, stringsAsFactors=FALSE)
            })
        })
        response$field_errors = ldply(response$objects, function(y) {
            ldply(y$status$field_errors, function(z) {
                z$resource = y$resource
                z$example = z$sample[[1]][[1]]
                z$count = z$sample[[1]][[2]]
                z$sample = NULL
                as.data.frame(z, stringsAsFactors=FALSE)
            })
        })
        response$row_format_errors = ldply(response$objects,
            function(y) {
                ldply(y$status$row_format_errors, function(z) {
                  z$resource = y$resource
                  x$status$field_errors = NULL
                  x$status$row_format_errors = NULL
                  z$example = z$sample[[1]][[1]]
                  z$count = z$sample[[1]][[2]]
                  z$sample = NULL
                  as.data.frame(z, stringsAsFactors=FALSE)
                })
            })
        response$objects = NULL
        if (datasets_only) {
            response = response$datasets
        }
    }
    return(response)
}
