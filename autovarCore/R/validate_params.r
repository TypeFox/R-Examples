#' Validates the params given to the autovar function
#'
#' This function uses a list of default params that may be overwritten the \code{params} argument. \code{stop()} errors are thrown when invalid params are supplied.
#' @param data_matrix The raw, unimputed data matrix. This parameter is supplied so that we can verify the selected column names.
#' @param params A \code{list} with the following named entries: \itemize{
#' \item \code{selected_column_names} - The endogenous variables in the models, specified as a vector of character strings. This argument is required. The selected column names should be a subset of the column names of \code{data_matrix}.
#' \item \code{significance_levels} - A vector with descending p values that indicate cut-offs placing models in different buckets. If it is not specified, this parameter defaults to \code{c(0.05, 0.01, 0.005)}. For example, with the default configuration, a model whose worst (lowest) p-level for any test is 0.03 is always seen as a better model than one whose worst p-level for any test is 0.009, no matter the AIC/BIC score of that model. Also, the lowest significance level indicates the minimum p-level for any test of a valid model. Thus, if a test for a model has a lower p-level than the minimum specified significance level, it is considered invalid.
#' \item \code{test_names} - The residual tests that should be performed, specified as a vector of character strings. If not specified, this parameter defaults to \code{c('portmanteau', 'portmanteau_squared', 'skewness')}. The possible tests are returned by the function \code{autovarCore:::supported_test_names()}. In addition to the residual tests, please note that the Eigenvalue stability test is always performed.
#' \item \code{criterion} - The information criterion used to sort the models. Valid options are 'AIC' (the default) or 'BIC'.
#' \item \code{imputation_iterations} - The number of times we average over one Amelia call for imputing the data set. Since one Amelia call averages over five imputations on its own, the actual number of imputations is five times the number specified here. The default value for this parameter is \code{30}.
#' \item \code{measurements_per_day} - The number of measurements per day in the time series data. The default value for this parameter is \code{1}.
#' }
#' @return A list containing augmented params.
#' @examples
#' data_matrix <- matrix(ncol = 3, nrow = 5)
#' data_matrix[, 1] <- 1
#' data_matrix[, 2] <- c(1, 3, 5, 6, 7)
#' data_matrix[, 3] <- c(1, 0, 1, NA, 1)
#' colnames(data_matrix) <- c('id', 'tijdstip', 'home')
#' autovarCore:::validate_params(data_matrix,
#'                               list(selected_column_names = c('tijdstip', 'home'),
#'                                    imputation_iterations = 20))
#' @export
validate_params <- function(data_matrix, params) {
  # precondition: dat_matrix is assumed to be a valid data set
  #               and is not validated here (it is validated in validate_raw_dataframe)
  returned_params <- default_autovar_params()
  assert_param_class(params, 'list')
  assert_param_subset(names(params),
                      c(names(returned_params), 'selected_column_names'))
  assert_param_presence('selected_column_names', names(params))
  for (param_name in names(params)) {
    validation_function <- switch(param_name,
       selected_column_names = validate_selected_column_names,
       significance_levels = validate_significance_levels,
       test_names = validate_test_names,
       criterion = validate_criterion,
       imputation_iterations = validate_imputation_iterations,
       measurements_per_day = validate_measurements_per_day)
    returned_params[[param_name]] <- validation_function(data_matrix,
                                                         params[[param_name]])
  }
  returned_params
}


# Validation functions

validate_selected_column_names <- function(data_matrix, given_param) {
  # precondition: data_matrix is a matrix
  assert_param_not_null(given_param)
  accepted_column_names <- colnames(data_matrix)
  assert_param_subset(given_param,
                      accepted_column_names,
                      "Invalid selected column name:")
  if (length(given_param) < 2)
    stop("Need at least two selected column names")
  if (length(given_param) > 31)
    stop("Need at most 31 selected column names")
  given_param
}

validate_significance_levels <- function(data_matrix, given_param) {
  assert_param_not_null(given_param)
  assert_param_class(given_param, 'numeric')
  sort(given_param, decreasing = TRUE)
}

validate_test_names <- function(data_matrix, given_param) {
  if (is.null(given_param)) return(NULL) # An empty vector of tests is allowed
  assert_param_subset(given_param,
                      supported_test_names(),
                      "Unsupported test name:")
  given_param
}

validate_criterion <- function(data_matrix, given_param) {
  assert_param_not_null(given_param)
  assert_param_single(given_param)
  assert_param_subset(given_param,
                      supported_criteria(),
                      "Unsupported criterion:")
  given_param
}

validate_imputation_iterations <- function(data_matrix, given_param) {
  assert_param_not_null(given_param)
  assert_param_single(given_param)
  assert_param_integer(given_param)
  assert_param_range(given_param, 1, 500, "number of imputation iterations")
  given_param
}

validate_measurements_per_day <- function(data_matrix, given_param) {
  assert_param_not_null(given_param)
  assert_param_single(given_param)
  assert_param_integer(given_param)
  # We may use the 0 value to denote that day- and daypart dummies are not to be included.
  assert_param_range(given_param, 0, 16, "number of measurements per day")
  given_param
}
