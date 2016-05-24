#' Creating BigML Models
#' @export
#' @family model methods
#' @references \url{https://bigml.com/developers/models}
#' @param dataset_id the relevant dataset_id used to create the model.
#' @param input_field_ids a vector of field ids to use for training.
#' @param name the name to give to the model.
#' @param objective_field_ids a vector of objective fields used for training.
#' @param range a vector of two values that define a range of instances from the dataset to train on.
#' @template dots
#' @return model_return
#' @template normal_methods
#' @examples
#' \dontrun{
#' # simple example
#' m1 = createModel("dataset/1")
#' # configure a number of different parameters
#' m2 = createModel("dataset/2", input_field_ids=c('000001'),
#'	objective_field_ids='000003', name='test', range = c(10,1000))
#' }
#' @references
#' \url{https://bigml.com/developers/datasets}
#' @template author

createModel <-
function (dataset_id, input_field_ids = NULL, name = NULL,
          objective_field_ids = NULL, range = NULL, ...)
{
    option = list()
    option$dataset = dataset_id
    if (!is.null(input_field_ids))
        option$input_fields = input_field_ids
    if (!is.null(name))
        option$name = name
    if (!is.null(objective_field_ids))
        option$objective_fields = objective_field_ids
    if (!is.null(range))
        option$range = range
    message("Model creation in progress...")
    return(.basic_api(.MODEL_URL)$postJson(option, ...))
}
