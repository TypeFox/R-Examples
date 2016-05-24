context('validate_params')

testdata_data_matrix <- function() {
  data_matrix <- matrix(ncol = 3, nrow = 5)
  data_matrix[, 1] <- 1
  data_matrix[, 2] <- c(1, 3, 5, 6, 7)
  data_matrix[, 3] <- c(1, 0, 1, NA, 1)
  colnames(data_matrix) <- c('id', 'tijdstip', 'home')
  data_matrix
}

test_that('validate_params requires the selected_column_names parameter', {
  expected_error_message <- "selected_column_names is a required parameter"
  data_matrix <- testdata_data_matrix()
  expect_error(autovarCore:::validate_params(data_matrix,
                                             list(measurements_per_day = 1)),
               expected_error_message)
  expect_error(autovarCore:::validate_params(data_matrix,
                                             list()),
               expected_error_message)
})

test_that('validate_params does not accept NULL for a params list', {
  expect_error(autovarCore:::validate_params(testdata_data_matrix(),
                                             NULL),
               "Param class should be: list")
})

test_that('validate_params substitutes with default parameter values', {
  expected_result <- autovarCore:::default_autovar_params()
  selected_column_names <- c('tijdstip', 'home')
  expected_result$selected_column_names <- selected_column_names
  expect_equal(autovarCore:::validate_params(testdata_data_matrix(),
                                             list(selected_column_names = selected_column_names)),
               expected_result)
})

test_that('validate_params does not accept invalid parameters', {
  expect_error(autovarCore:::validate_params(testdata_data_matrix(),
                                             list(selected_column_names = c('tijdstip', 'home'),
                                                  invalid_param = 'something')),
               "Invalid param: invalid_param")
})


test_that('validate_params calls the correct subfunctions to override the default parameters', {
  expected_result <- list(significance_levels = 1,
                          test_names = 2,
                          criterion = 3,
                          imputation_iterations = 4,
                          measurements_per_day = 5,
                          selected_column_names = 6)
  with_mock(
    `autovarCore:::validate_significance_levels` = function(...) 1,
    `autovarCore:::validate_test_names` = function(...) 2,
    `autovarCore:::validate_criterion` = function(...) 3,
    `autovarCore:::validate_imputation_iterations` = function(...) 4,
    `autovarCore:::validate_measurements_per_day` = function(...) 5,
    `autovarCore:::validate_selected_column_names` = function(...) 6,
     expect_equal(autovarCore:::validate_params(testdata_data_matrix(), expected_result),
                  expected_result)
  )
})


# Validation functions

test_that('validate_selected_column_names accepts only names of columns in the data frame', {
  data_matrix <- testdata_data_matrix()
  expect_error(autovarCore:::validate_selected_column_names(data_matrix, NULL),
               "Given param cannot be NULL")
  expect_error(autovarCore:::validate_selected_column_names(data_matrix,
                                                            'unknown'),
               "Invalid selected column name: unknown")

  expect_error(autovarCore:::validate_selected_column_names(data_matrix,
                                                            c('tijdstip', 'unknown')),
               "Invalid selected column name: unknown")
  # Below goes wrong because the minimum required number of columns is 2.
  expect_error(autovarCore:::validate_selected_column_names(data_matrix,
                                                            'tijdstip'),
               "Need at least two selected column names")
  col_names <- c('aa', 'ab', 'ac', 'ad', 'ae', 'af', 'ag', 'ah',
                 'ai', 'aj', 'ak', 'al', 'am', 'an', 'ao', 'ap',
                 'aq', 'ar', 'as', 'at', 'au', 'av', 'aw', 'ax',
                 'ay', 'az', 'ba', 'bb', 'bc', 'bd', 'be', 'bf')
  data_matrix2 <- matrix(NA, nrow = 1, ncol = 32, dimnames = list(NULL, col_names))
  expect_error(autovarCore:::validate_selected_column_names(data_matrix2,
                                                            col_names),
               "Need at most 31 selected column names")
  # The statements below should not throw errors.
  expect_equal(autovarCore:::validate_selected_column_names(data_matrix,
                                                            c('tijdstip', 'home')),
               c('tijdstip', 'home'))
  expect_equal(autovarCore:::validate_selected_column_names(data_matrix,
                                                            c('tijdstip', 'home', 'id')),
               c('tijdstip', 'home', 'id'))
})

test_that('validate_significance_levels returns the input sorted decreasingly', {
  expect_equal(autovarCore:::validate_significance_levels(testdata_data_matrix(),
                                                          c(0.5, 0.3, 0.7)),
                c(0.7, 0.5, 0.3))
})

test_that('validate_significance_levels accepts only numeric vectors', {
  data_matrix <- testdata_data_matrix()
  expect_error(autovarCore:::validate_significance_levels(data_matrix, NULL),
               "Given param cannot be NULL")
  expect_error(autovarCore:::validate_significance_levels(data_matrix, c('a', 'b', 'c')),
               "Param class should be: numeric")
  # The statement below should not throw an error.
  expect_equal(autovarCore:::validate_significance_levels(data_matrix, 0.7),
                0.7)
})

test_that('validate_test_names accepts only supported test names', {
  data_matrix <- testdata_data_matrix()
  expect_error(autovarCore:::validate_test_names(data_matrix, c('portmanteau', 'unknown')),
               "Unsupported test name: unknown")
  expect_error(autovarCore:::validate_test_names(data_matrix, 'unknown'),
               "Unsupported test name: unknown")
  # The statements below should not throw errors.
  expect_null(autovarCore:::validate_test_names(data_matrix, NULL))
  expect_equal(autovarCore:::validate_test_names(data_matrix, c('portmanteau', 'skewness')),
               c('portmanteau', 'skewness'))
  expect_equal(autovarCore:::validate_test_names(data_matrix, 'kurtosis'),
               'kurtosis')
})

test_that('validate_criterion accepts only supported criteria', {
  data_matrix <- testdata_data_matrix()
  expect_error(autovarCore:::validate_criterion(data_matrix, NULL),
               "Given param cannot be NULL")
  expect_error(autovarCore:::validate_criterion(data_matrix, 'unknown'),
               "Unsupported criterion: unknown")
  expect_error(autovarCore:::validate_criterion(data_matrix, c('AIC', 'BIC')),
               "Length of given param is not 1:")
  # The statement below should not throw an error.
  expect_equal(autovarCore:::validate_criterion(data_matrix, 'AIC'),
               'AIC')
  expect_equal(autovarCore:::validate_criterion(data_matrix, 'BIC'),
               'BIC')
})

test_that('validate_imputation_iterations accepts only an integer in range', {
  data_matrix <- testdata_data_matrix()
  expect_error(autovarCore:::validate_imputation_iterations(data_matrix, NULL),
               "Given param cannot be NULL")
  expect_error(autovarCore:::validate_imputation_iterations(data_matrix, c(2, 3)),
               "Length of given param is not 1:")
  expect_error(autovarCore:::validate_imputation_iterations(data_matrix, list(a = 2, b = 3)),
               "Length of given param is not 1:")
  expect_error(autovarCore:::validate_imputation_iterations(data_matrix, 3.5),
               "Given param is not an integer:")
  expect_error(autovarCore:::validate_imputation_iterations(data_matrix, 'hoi'),
               "Given param is not an integer:")
  expect_error(autovarCore:::validate_imputation_iterations(data_matrix, 0),
               "The number of imputation iterations has to be an integer in range 1-500")
  expect_error(autovarCore:::validate_imputation_iterations(data_matrix, 501),
               "The number of imputation iterations has to be an integer in range 1-500")
  # The statements below should not throw errors.
  expect_equal(autovarCore:::validate_imputation_iterations(data_matrix, 1), 1)
  expect_equal(autovarCore:::validate_imputation_iterations(data_matrix, 500), 500)
  expect_equal(autovarCore:::validate_imputation_iterations(data_matrix, 376), 376)
})

test_that('validate_measurements_per_day accepts only an integer in range', {
  data_matrix <- testdata_data_matrix()
  expect_error(autovarCore:::validate_measurements_per_day(data_matrix, NULL),
               "Given param cannot be NULL")
  expect_error(autovarCore:::validate_measurements_per_day(data_matrix, c(2, 3)),
               "Length of given param is not 1:")
  expect_error(autovarCore:::validate_measurements_per_day(data_matrix, list(a = 2, b = 3)),
               "Length of given param is not 1:")
  expect_error(autovarCore:::validate_measurements_per_day(data_matrix, 3.5),
               "Given param is not an integer:")
  expect_error(autovarCore:::validate_measurements_per_day(data_matrix, 'hoi'),
               "Given param is not an integer:")
  expect_error(autovarCore:::validate_measurements_per_day(data_matrix, -1),
               "The number of measurements per day has to be an integer in range 0-16")
  expect_error(autovarCore:::validate_measurements_per_day(data_matrix, 17),
               "The number of measurements per day has to be an integer in range 0-16")
  # The statements below should not throw errors.
  expect_equal(autovarCore:::validate_measurements_per_day(data_matrix, 0), 0)
  expect_equal(autovarCore:::validate_measurements_per_day(data_matrix, 13), 13)
  expect_equal(autovarCore:::validate_measurements_per_day(data_matrix, 16), 16)
})
