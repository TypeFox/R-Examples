#' Calculate dummy variables to mask residual outliers
#'
#' This function returns a matrix with columns that have a \code{1} at indices where the residuals have an outlier, and a \code{0} everywhere else. Outliers are calculated per variable (column) separately. We consider residual outliers the rows in the column of residuals or in the column of squared residuals that are more than 2.5 times the standard deviation away from the mean (standard deviation and mean are calculated separately per column and for residuals/squared residuals). The dummy columns are prepended with zeros to match the size of the other input variables to the model.
#' @param resid_matrix A matrix of residuals. Column names are copied over to the returned result.
#' @param number_of_rows The number of measurements that were input to the model. Since the length of the residual matrix is shorter depending on the amount of lags in the model, we use \code{number_of_rows} to specify the number of rows in the returned matrix.
#' @return A matrix with dummy variables in columns following the procedure described above.
#' @examples
#' resid_matrix <- matrix(rnorm(39 * 3),
#'                        nrow = 39,
#'                        ncol = 3,
#'                        dimnames = list(NULL, c('rumination', 'happiness', 'activity')))
#' resid_matrix[13, 2] <- 48
#' resid_matrix[23, 2] <- -62
#' resid_matrix[36, 2] <- 33
#' resid_matrix[27, 3] <- 75
#' resid_matrix
#' autovarCore:::residual_outliers(resid_matrix, 40)
#' @export
residual_outliers <- function(resid_matrix, number_of_rows) {
  result <- matrix(NA,
                   ncol = ncol(resid_matrix),
                   nrow = number_of_rows,
                   dimnames = list(NULL, colnames(resid_matrix)))
  for (column_index in 1:ncol(resid_matrix)) {
    result[, column_index] <- residual_outliers_column(resid_matrix[, column_index],
                                                       number_of_rows)
  }
  result
}

residual_outliers_column <- function(resid_column, number_of_rows) {
  result <- rep.int(0, number_of_rows)
  result <- pmax(result, normal_outliers_column(resid_column, number_of_rows))
  result <- pmax(result, squared_outliers_column(resid_column, number_of_rows))
  result
}

normal_outliers_column <- function(resid_column, number_of_rows) {
  outliers_column(resid_column,
                  number_of_rows,
                  std_factor_for_normal_outliers())
}

squared_outliers_column <- function(resid_column, number_of_rows) {
  outliers_column(resid_column * resid_column,
                  number_of_rows,
                  std_factor_for_squared_outliers())
}

outliers_column <- function(column_data, number_of_rows, std_factor) {
  result <- as.numeric(abs(column_data - mean(column_data)) > std_factor * sd(column_data))
  if (length(result) < number_of_rows)
    result <- c(rep.int(0, number_of_rows - length(result)), result)
  result
}
