#' Quickly Creating BigML Predictions
#' @export
#' @family prediction methods
#' @references \url{https://bigml.com/developers/predictions}
#' @family quick methods
#' @param model A character string or response object containing a valid model
#' id value.
#' @param values A named vector or list of elements to retrieve a prediction
#'	for
#' @param name A string giving the name of the prediction.
#' @param prediction_only if TRUE, only the predicted value is returned.
#' 	Otherwise, the full API response is returned.
#' @template dots
#' @template prediction_return
#' @template author
#' @return A numeric or string value giving the prediction.
#' @details quickPrediction can operate on a model id string, or a model
#'	response object from an earlier request.  The \code{values} are a list of
#'	named elements that are used as input.
#' @examples
#' \dontrun{
#' quickPrediction("model/1", list(Sepal.Width=3.5, Petal.Length=1.4))
#' # 'setosa'
#' }
quickPrediction <-
function (model, values, name = NULL, prediction_only = TRUE, ...)
{
    model_id = .resolve_resource_id(model, "model")
    if (is.null(names(values))) {
        stop("values argument must have named values")
    }
    option = list()
    option$model = model_id
    mresponse = getModel(model_id, ...)
    idlist = list()
    for (v in names(values)) {
        m_id = .resolve_field_id(v, mresponse$model$fields)
        if (is.null(m_id)) {
            stop(paste("value does not exist in model:", v))
        }
        idlist[m_id] = values[v]
    }
    option$input_data = idlist
    if (!is.null(name))
        option$name = name
    response = .basic_api(.PREDICTION_URL)$postJson(option, ...)
    if (prediction_only)
        return(as.vector(response$prediction))
    else
        return(response)
}
