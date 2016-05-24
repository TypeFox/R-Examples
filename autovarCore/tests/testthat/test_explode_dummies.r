context('explode_dummies')

test_that('explode_dummies calls merge_columns correctly', {
  called_count_merge_columns <<- 0
  called_count_outlier_dummy_column <<- 0
  input_matrix <<- matrix(NA, nrow = 3, ncol = 2, dimnames = list(NULL, c('a', 'b')))
  input_matrix[, 1] <<- c(0, 0, 1)
  input_matrix[, 2] <<- c(0, 1, 0)
  expected_result <- matrix(1, nrow = 1, ncol = 3)
  with_mock(
    `autovarCore:::merge_columns` = function(...) {
      called_count_merge_columns <<- called_count_merge_columns + 1
      expect_equal(list(...), list(input_matrix))
      c(1, 1, 1)
    },
    `autovarCore:::outlier_dummy_column` = function(...) {
      called_count_outlier_dummy_column <<- called_count_outlier_dummy_column + 1
      expect_equal(list(...), list(called_count_outlier_dummy_column, 3))
      1
    },
    expect_equal(autovarCore:::explode_dummies(input_matrix), expected_result)
  )
  expect_equal(called_count_merge_columns, 1)
  expect_equal(called_count_outlier_dummy_column, 3)
  rm(list = c('called_count_merge_columns',
              'called_count_outlier_dummy_column',
              'input_matrix'), pos = '.GlobalEnv')
})

test_that('explode_dummies returns the correct result', {
  input_matrix <- matrix(NA, nrow = 6, ncol = 3, dimnames = list(NULL, c('a', 'b', 'c')))
  input_matrix[, 1] <- c(0, 0, 1, 0, 1, 1)
  input_matrix[, 2] <- c(1, 0, 1, 0, 0, 0)
  input_matrix[, 3] <- c(0, 1, 1, 0, 1, 0)
  expected_result <- matrix(NA, nrow = 6, ncol = 5,
                            dimnames = list(NULL, c('outlier_1',
                                                    'outlier_2',
                                                    'outlier_3',
                                                    'outlier_5',
                                                    'outlier_6')))
  expected_result[, 1] <- c(1, 0, 0, 0, 0, 0)
  expected_result[, 2] <- c(0, 1, 0, 0, 0, 0)
  expected_result[, 3] <- c(0, 0, 1, 0, 0, 0)
  expected_result[, 4] <- c(0, 0, 0, 0, 1, 0)
  expected_result[, 5] <- c(0, 0, 0, 0, 0, 1)
  expect_equal(autovarCore:::explode_dummies(input_matrix), expected_result)
})


test_that('merge_columns returns the correct result', {
  input_matrix <- matrix(NA, nrow = 5, ncol = 3, dimnames = list(NULL, c('a', 'b', 'c')))
  input_matrix[, 1] <- c(0, 0, 1, 0, 1)
  input_matrix[, 2] <- c(1, 0, 1, 0, 0)
  input_matrix[, 3] <- c(0, 1, 1, 0, 1)
  expected_result <- matrix(NA, nrow = 5, ncol = 4,
                            dimnames = list(NULL, c('outlier_1',
                                                    'outlier_2',
                                                    'outlier_3',
                                                    'outlier_5')))
  expected_result <- c(1, 1, 1, 0, 1)
  expect_equal(autovarCore:::merge_columns(input_matrix), expected_result)
})

test_that('outlier_dummy_column returns the correct result', {
  expected_result <- matrix(NA, nrow = 5, ncol = 1, dimnames = list(NULL, 'outlier_1'))
  expected_result[, 1] <- c(1, 0, 0, 0, 0)
  expect_equal(autovarCore:::outlier_dummy_column(1, 5), expected_result)
  expected_result <- matrix(NA, nrow = 6, ncol = 1, dimnames = list(NULL, 'outlier_6'))
  expected_result[, 1] <- c(0, 0, 0, 0, 0, 1)
  expect_equal(autovarCore:::outlier_dummy_column(6, 6), expected_result)
  expected_result <- matrix(NA, nrow = 4, ncol = 1, dimnames = list(NULL, 'outlier_2'))
  expected_result[, 1] <- c(0, 1, 0, 0)
  expect_equal(autovarCore:::outlier_dummy_column(2, 4), expected_result)
})
