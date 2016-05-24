context('residual_outliers')

test_that('residual_outliers returns the correct result', {
  var_names <- c('a', 'b', 'c')
  expected_result <- matrix(0, ncol = 3, nrow = 10,
                            dimnames = list(NULL, var_names))
  expected_result[2, 1] <- 1
  expected_result[7, 2] <- 1
  expected_result[10, 3] <- 1
  input_matrix <- matrix(NA, ncol = 3, nrow = 9,
                         dimnames = list(NULL, var_names))
  input_matrix[, 1] <- c(20, 0.2, 0.4, 0.3, -0.3, 0.2, 0.1, 0.3, 0.5)
  input_matrix[, 2] <- c(0.2, 0.2, 0.4, 0.3, -0.3, 20, 0.1, 0.3, 0.5)
  input_matrix[, 3] <- c(-0.1, 0.2, 0.4, 0.2, -0.3, 0.2, 0.1, 0.3, -50)
  expect_equal(autovarCore:::residual_outliers(input_matrix, 10),
               expected_result)
})

test_that('residual_outliers works with different values for number_of_rows', {
  var_names <- c('a', 'b', 'c')
  expected_result <- matrix(0, ncol = 3, nrow = 15,
                            dimnames = list(NULL, var_names))
  expected_result[7, 1] <- 1
  expected_result[12, 2] <- 1
  expected_result[15, 3] <- 1
  input_matrix <- matrix(NA, ncol = 3, nrow = 9,
                         dimnames = list(NULL, var_names))
  input_matrix[, 1] <- c(20, 0.2, 0.4, 0.3, -0.3, 0.2, 0.1, 0.3, 0.5)
  input_matrix[, 2] <- c(0.2, 0.2, 0.4, 0.3, -0.3, 20, 0.1, 0.3, 0.5)
  input_matrix[, 3] <- c(-0.1, 0.2, 0.4, 0.2, -0.3, 0.2, 0.1, 0.3, -50)
  expect_equal(autovarCore:::residual_outliers(input_matrix, 15),
               expected_result)
})

test_that('residual_outliers calls residual_outlier_column correctly', {
  called_count <<- 0
  var_names <- c('a', 'b')
  input_matrix <- matrix(NA, ncol = 2, nrow = 2, dimnames = list(NULL, var_names))
  expected_result <- matrix(1, ncol = 2, nrow = 5, dimnames = list(NULL, var_names))
  with_mock(
    `autovarCore:::residual_outliers_column` = function(...) {
      called_count <<- called_count + 1
      expect_equal(list(...), list(rep(NA, 2), 5))
      rep.int(1, 5)
    },
    expect_equal(autovarCore:::residual_outliers(input_matrix, 5),
                 expected_result)
  )
  expect_equal(called_count, 2)
  rm(list = 'called_count', pos = '.GlobalEnv')
})


test_that('residual_outliers_column combines the results of its dependent functions correctly', {
  called_count_normal <<- 0
  called_count_squared <<- 0
  with_mock(
    `autovarCore:::normal_outliers_column` = function(...) {
      called_count_normal <<- called_count_normal + 1
      expect_equal(list(...), list(1, 2))
      c(1, 0)
    },
    `autovarCore:::squared_outliers_column` = function(...) {
      called_count_squared <<- called_count_squared + 1
      expect_equal(list(...), list(1, 2))
      c(0, 1)
    },
    expect_equal(autovarCore:::residual_outliers_column(1, 2), c(1, 1))
  )
  expect_equal(called_count_normal, 1)
  expect_equal(called_count_squared, 1)
  rm(list = c('called_count_normal',
              'called_count_squared'), pos = '.GlobalEnv')
})


test_that('normal_outliers_column calls outliers_column and std_factor_for_normal_outliers correclty', {
  called_count_std_factor <<- 0
  called_count_outliers_column <<- 0
  with_mock(
    `autovarCore:::std_factor_for_normal_outliers` = function(...) {
      called_count_std_factor <<- called_count_std_factor + 1
      expect_equal(list(...), list())
      3
    },
    `autovarCore:::outliers_column` = function(...) {
      called_count_outliers_column <<- called_count_outliers_column + 1
      expect_equal(list(...), list(1, 2, 3))
      4
    },
    expect_equal(autovarCore:::normal_outliers_column(1, 2), 4)
  )
  expect_equal(called_count_std_factor, 1)
  expect_equal(called_count_outliers_column, 1)
  rm(list = c('called_count_std_factor',
              'called_count_outliers_column'), pos = '.GlobalEnv')
})

test_that('squared_outliers_column calls outliers_column and std_factor_for_squared_outliers correclty', {
  called_count_std_factor <<- 0
  called_count_outliers_column <<- 0
  with_mock(
    `autovarCore:::std_factor_for_squared_outliers` = function(...) {
      called_count_std_factor <<- called_count_std_factor + 1
      expect_equal(list(...), list())
      30
    },
    `autovarCore:::outliers_column` = function(...) {
      called_count_outliers_column <<- called_count_outliers_column + 1
      expect_equal(list(...), list(10 * 10, 20, 30))
      40
    },
    expect_equal(autovarCore:::squared_outliers_column(10, 20), 40)
  )
  expect_equal(called_count_std_factor, 1)
  expect_equal(called_count_outliers_column, 1)
  rm(list = c('called_count_std_factor',
              'called_count_outliers_column'), pos = '.GlobalEnv')
})


test_that('outlier_column calculates the outliers correctly', {
  column_data <- c(0.5, 0.2, 0.4, 0.3, -5.2, 0.2, 0.1, 0.3, 20)
  expected_result <- c(0, 0, 0, 0, 0, 0, 0, 0, 1)
  expect_equal(autovarCore:::outliers_column(column_data,
                                             length(column_data),
                                             2.5),
               expected_result)
  expected_result <- c(0, 0, 0, 0, 1, 0, 0, 0, 1)
  expect_equal(autovarCore:::outliers_column(column_data,
                                             length(column_data),
                                             1),
               expected_result)

})

test_that('outlier_column respects values for number_of_rows', {
  column_data <- c(0.2, -0.1, 0.4, -0.4, 33, 0.3, -0.2, -0.4, 0.2, -0.3, 0.4)
  expected_result <- c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0)
  expect_equal(autovarCore:::outliers_column(column_data, 11, 2.5), expected_result)
  expected_result <- c(0, expected_result)
  expect_equal(autovarCore:::outliers_column(column_data, 12, 2.5), expected_result)
  expected_result <- c(rep.int(0, 5), expected_result)
  expect_equal(autovarCore:::outliers_column(column_data, 17, 2.5), expected_result)
})
