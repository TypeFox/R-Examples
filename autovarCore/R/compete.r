#' Returns the winning model
#'
#' This function returns the best model as explained in the documentation for the \code{autovar} function.
#' @param best A model given as a list with at least the properties \code{model_score, nr_dummy_variables,} and \code{bucket}.
#' @param challenger Another model, also given as a list with properties \code{model_score, nr_dummy_variables,} and \code{bucket}.
#' @param compare_outliers A boolean. When \code{FALSE}, the model comparison does not take the number of dummy variables into account.
#' @return This function returns the best model of the two given models.
#' @examples
#' model1 <- list(logtransformed = FALSE, lag = 1, nr_dummy_variables = 1,
#'                model_score = 100, bucket = 0.05)
#' model2 <- list(logtransformed = FALSE, lag = 2, nr_dummy_variables = 2,
#'                model_score = 200, bucket = 0.01)
#' autovarCore:::compete(model1, model2, TRUE)
#' @export
compete <- function(best, challenger, compare_outliers) {
  if (challenger_wins(best, challenger, compare_outliers))
    return(challenger)
  best
}

challenger_wins <- function(best, challenger, compare_outliers) {
  if (challenger$bucket != best$bucket)
    return(challenger$bucket > best$bucket)
  if (compare_outliers && challenger$nr_dummy_variables != best$nr_dummy_variables)
    return(challenger$nr_dummy_variables < best$nr_dummy_variables)
  challenger$model_score < best$model_score
}
