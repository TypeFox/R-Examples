#' Calculate a bit mask to identify invalid outlier dummies
#'
#' Invalid outlier dummy variables are dummy variables that are all zeros (where the original variable had no outliers at the 2.5x standard deviation for either the residuals or the squared residuals. Interpreting the leftmost column as bit 0 and continuing with higher bits going from left to right in the matrix, this function returns a bit mask that has a 1 on all positions in the matrix where the dummy column is invalid. We use this in later functions to easily filter out theinvalid outlier masks from the valid ones.
#' @param outlier_dummies A matrix of outlier dummy variables in columns.
#' @return An integer mask indicating the invalid columns according to the procedure describe above.
#' @examples
#' resid_matrix <- matrix(rnorm(39 * 3),
#'                        nrow = 39,
#'                        ncol = 3,
#'                        dimnames = list(NULL, c('rumination', 'happiness', 'activity')))
#' outlier_dummies <- autovarCore:::residual_outliers(resid_matrix, 40)
#' autovarCore:::invalid_mask(outlier_dummies)
#' @export
invalid_mask <- function(outlier_dummies) {
  result <- 0
  for (column_index in 1:ncol(outlier_dummies))
    if (column_is_invalid(outlier_dummies[, column_index]))
      result <- result + bitwShiftL(1, column_index - 1)
  result
}

column_is_invalid <- function(outlier_column) {
  all(outlier_column == 0)
}
