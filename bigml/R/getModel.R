#' Retrieving a BigML Model
#' @export
#' @family model methods
#' @references \url{https://bigml.com/developers/models}
#' @param model_id A string giving the model id.
#' @template dots
#' @template model_return
#' @template normal_methods
#' @template author
getModel <-
function (model_id, ...)
{
    message("Retrieving the model...")
    return(.basic_api(.MODEL_URL)$get(id=model_id, ...))
}
