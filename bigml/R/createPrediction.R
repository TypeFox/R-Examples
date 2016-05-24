#' Creating BigML Predictions
#' @export
#' @family prediction methods
#' @references \url{https://bigml.com/developers/predictions}
#' @param model_id character string; the model id
#' @param input_field_ids a list of input field ids and values to make a
#'	prediction for (see example).
#' @param name character string; The given name for the prediction.
#' @param prediction_only logical: Indicating whether the prediction should
#' be returned as a simple value, or if the full response object should be
#' returned.
#' @template dots
#' @template prediction_return
#' @template normal_methods
#' @examples
#' \dontrun{
#' # simple example
#' m1 = createPrediction("model/1",
#'	input_field_ids = c('000001'='somevalue', '000002'=9999))
#' # configure a number of different parameters
#' m2 = createPrediction("model/2",
#'	input_field_ids = c('000001'='somevalue', '000002'=9999),
#'	name='new prediction')
#' }
#' @template author
createPrediction <-
function (model_id, input_field_ids, name = NULL, prediction_only=TRUE, ...)
{
    option = list()
    option$input_data = input_field_ids
	option$model = model_id
    if (!is.null(name))
        option$name = name
    message("Prediction creation in progress...")
    response = .basic_api(.PREDICTION_URL)$post(option, ...)
    if (prediction_only)
        return(as.vector(response$prediction))
    else
	    return(response)
}
