#' Listing BigML Datasets
#' @export
#' @family model methods
#' @references \url{https://bigml.com/developers/models}
#' @param flatten A logical value indicating whether to flatten the response
#'	into a data frame.
#' @param models_only A logical value indicating whether to only return
#'	the data frame of model information (only valid if \code{flatten} is
#'	\code{TRUE}).
#' @template dots
#' @template model_list_return
#' @template author
listModels <-
function (flatten = TRUE, models_only = TRUE, ...)
{
    message("Retrieving the model list...")
    response = .basic_api(.MODEL_URL)$list(...)
    if (flatten) {
        response$models = ldply(response$objects, function(x) {
            as.data.frame(Filter(function(y) length(y) == 1,
                x), stringsAsFactors=FALSE)
        })
		response$objects = NULL
    }
    if (models_only) {
        response = response$models
    }
    response
}
