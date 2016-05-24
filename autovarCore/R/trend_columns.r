#' Construct linear and quadratic trend columns
#'
#' This function returns a matrix of linear and quadratic trends.
#' @param number_of_rows the number of rows in the input data set.
#' @return A matrix with \code{number_of_rows} rows and 2 columns, one for linear trends and one for quadratic trends.
#' @examples
#' autovarCore:::trend_columns(10)
#' @export
trend_columns <- function(number_of_rows) {
  result <- NULL
  result <- add_linear_trend(result, number_of_rows)
  result <- add_squared_trend(result, number_of_rows)
  result
}

add_linear_trend <- function(column_matrix, number_of_rows) {
  column <- matrix(1:number_of_rows, dimnames = list(NULL, 'index'))
  cbind(column_matrix, column)
}

add_squared_trend <- function(column_matrix, number_of_rows) {
  column <- matrix((1:number_of_rows)^2, dimnames = list(NULL, 'index2'))
  cbind(column_matrix, column)
}
