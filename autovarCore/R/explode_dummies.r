#' Explode dummies columns into separate dummy variables
#'
#' This function takes a matrix with dummy outlier columns, where there are possibly multiple ones. We first merge these columns to one and then explode them to obtain one dummy variable per column.
#' @param outlier_dummies A matrix of outlier dummy variables in columns.
#' @return A matrix with dummy variables in columns, each having one nonzero index. The columns are named \code{outlier_x}, with x being the 1-based row index of the position that this dummy variable is masking.
#' @examples
#' outlier_dummies <- matrix(NA,
#'                           nrow = 5,
#'                           ncol = 3,
#'                          dimnames = list(NULL, c('rumination', 'happiness', 'activity')))
#' outlier_dummies[, 1] <- c(0, 0, 1, 0, 1)
#' outlier_dummies[, 2] <- c(0, 1, 1, 0, 0)
#' outlier_dummies[, 3] <- c(1, 0, 0, 0, 1)
#' outlier_dummies
#' autovarCore:::explode_dummies(outlier_dummies)
#' @export
explode_dummies <- function(outlier_dummies) {
  merged_vector <- merge_columns(outlier_dummies)
  result <- NULL
  number_of_rows <- nrow(outlier_dummies)
  for (row_index in 1:number_of_rows)
    if (merged_vector[row_index] == 1)
      result <- cbind(result, outlier_dummy_column(row_index, number_of_rows))
  result
}

merge_columns <- function(outlier_dummies) {
  result <- NULL
  for (row_index in 1:nrow(outlier_dummies))
    result <- c(result, max(outlier_dummies[row_index, ]))
  result
}

outlier_dummy_column <- function(row_index, number_of_rows) {
  result <- matrix(0, ncol = 1, nrow = number_of_rows,
                   dimnames = list(NULL, paste('outlier_', row_index, sep = '')))
  result[row_index, ] <- 1
  result
}
