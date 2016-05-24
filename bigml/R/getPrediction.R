#' Retrieving a BigML Prediction
#' @export
#' @family prediction methods
#' @references \url{https://bigml.com/developers/predictions}
#' @param prediction_id the id of the prediction resource.
#' @template dots
#' @template prediction_return
#' @template normal_methods
#' @template author
getPrediction <-
function (prediction_id, ...)
{
    message("Retrieving the prediction...")
    return(.basic_api(.PREDICTION_URL)$get(id = prediction_id))
}