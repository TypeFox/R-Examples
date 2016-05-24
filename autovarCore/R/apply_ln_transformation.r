#' Applies the natural logarithm to the data set
#'
#' This applies the ln function columnwise to the given input matrix and returns the modified matrix. If necessary, columns undergo a linear translation to ensure that all resulting values are >= 0.
#' @param data_matrix The original data matrix.
#' @return The log-transformed data matrix.
#' @examples
#' data_matrix <- matrix(1:10, dimnames = list(NULL, 'some_val'))
#' data_matrix
#' autovarCore:::apply_ln_transformation(data_matrix)
#' @export
apply_ln_transformation <- function(data_matrix) {
  for (column_index in 1:ncol(data_matrix))
    data_matrix[, column_index] <- ln_column(data_matrix[, column_index])
  data_matrix
}
ln_column <- function(data_column) {
  # Note that in the current use of apply_ln_transformation,
  # the data_matrix will not contain any NA values. The na.rm
  # argument below is just so this function will work in a
  # broader context later on, if needed.
  increment <- 1 - min(data_column, na.rm = TRUE)
  if (increment > 0)
    data_column <- data_column + increment
  log(data_column)
}
