#' Imputes the missing values in the input data
#'
#' This function uses \code{Amelia::amelia} to try and impute missing (\code{NA}) values in the input data set. Amelia averages over five iterations, which is multiplied by the given \code{imputation_iterations} parameter.
#' @param data_matrix The raw, unimputed data matrix.
#' @param measurements_per_day The number of measurements per day. This variable is used for adding day part dummy variables to aid the imputation.
#' @param imputation_iterations The amount of times the Amelia call should be averaged over. The actual number of imputations is five times the value for \code{imputation_iterations}, since Amelia's values are already averaged over five runs.
#' @return This function returns the modified matrix.
#' @examples
#' # create a matrix with some missing values
#' data_matrix <- matrix(nrow = 40, ncol = 3)
#' data_matrix[, ] <- runif(ncol(data_matrix) * nrow(data_matrix), 1, nrow(data_matrix))
#' while (sum(is.na(data_matrix)) == 0)
#'   data_matrix[as.logical(round(runif(ncol(data_matrix) * nrow(data_matrix), -0.3, 0.7)))] <- NA
#' colnames(data_matrix) <- c('rumination', 'happiness', 'activity')
#' data_matrix
#' autovarCore:::impute_datamatrix(data_matrix, 1, 30)
#' @export
impute_datamatrix <- function(data_matrix, measurements_per_day, imputation_iterations) {
  # precondition: imputation_iterations > 0. This precondition is met by validate_params.
  if (!has_missings(data_matrix))
    return(data_matrix)
  orig_number_of_columns <- ncol(data_matrix)
  data_matrix <- cbind(data_matrix,
                       time = 1:nrow(data_matrix),
                       daypart = rep(0:(measurements_per_day - 1),
                                     nrow(data_matrix))[1:nrow(data_matrix)])
  constant_columns <- NULL
  if (length(colnames(data_matrix)) != dim(data_matrix)[[2]])
    stop("Unnamed columns found in matrix")
  for (column_name in colnames(data_matrix))
    if (all(is.na(data_matrix[, column_name])) ||
        is.na(var(data_matrix[, column_name], na.rm = TRUE)) ||
        var(data_matrix[, column_name], na.rm = TRUE) == 0) {
      constant_columns <- c(constant_columns, column_name)
      # If there is no variance but some rows still have a
      # value that is not NA, set the whole column to that value.
      if (!is.na(mean(data_matrix[, column_name], na.rm = TRUE)))
        data_matrix[, column_name] <- mean(data_matrix[, column_name], na.rm = TRUE)
    }
  nominal_variables <- 'daypart'
  if ('daypart' %in% constant_columns)
    nominal_variables <- NULL
  variable_part <- as.matrix(data_matrix[, !(colnames(data_matrix) %in% constant_columns)])
  tolerance <- 0.1
  time_variable_name <- 'time'
  column_numbers_of_lagged_variables <- 1:(ncol(variable_part) - 1 -
                                        ifelse(is.null(nominal_variables), 0, 1))
  power_of_the_imputation_polynomial <- 2
  output_detail <- 0 # no screen output. 1 for normal output, 2 for verbose output
  level_of_empirical_prior <- 0.01 * nrow(variable_part)
  allow_increasing_empirical_prior <- 1
  imputed_variable_part <- 0
  successful_imputations <- 0
  for (i in 1:imputation_iterations) {
    amelia_result <- amelia(variable_part,
                            tol = tolerance,
                            ts = time_variable_name,
                            lags = column_numbers_of_lagged_variables,
                            noms = nominal_variables,
                            polytime = power_of_the_imputation_polynomial,
                            p2s = output_detail,
                            empri = level_of_empirical_prior,
                            autopri = allow_increasing_empirical_prior)
    if (is.null(amelia_result$imputations))
      next
    successful_imputations <- successful_imputations + 1
    mean_imputation_result <- .2 * Reduce('+', amelia_result$imputations)
    imputed_variable_part <- imputed_variable_part + mean_imputation_result
  }
  if (successful_imputations == 0)
    stop("Amelia imputation failed")
  imputed_variable_part = imputed_variable_part / successful_imputations
  data_matrix[, !(colnames(data_matrix) %in% constant_columns)] <- imputed_variable_part
  data_matrix[, 1:orig_number_of_columns]
}

has_missings <- function(data_matrix) {
  sum(is.na(data_matrix)) > 0
}
