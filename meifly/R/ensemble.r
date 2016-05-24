#' Create a new ensemble of models.
#'
#' @keywords internal
#' @param models list of models
#' @param data associated data frame
new_ensemble <- function(models, data) {
  structure(models, data = data, class = "ensemble")
}

#' @export
"[.ensemble" <- function(x, i, ...) {
  structure(NextMethod(), class = "ensemble")
}
