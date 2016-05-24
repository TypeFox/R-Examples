context('invalid_mask')

test_that('invalid_mask returns the correct result', {
  input_matrix <- matrix(NA, nrow = 5, ncol = 3,
                          dimnames = list(NULL, c('rumination', 'happiness', 'activity')))
  input_matrix[, 'rumination'] <- c(1, 0, 0, 0, 0)
  input_matrix[, 'happiness'] <- c(0, 1, 0, 0, 0)
  input_matrix[, 'activity'] <- c(0, 0, 0, 0, 0)
  expect_equal(autovarCore:::invalid_mask(input_matrix), 4)
})

test_that('invalid_mask adds up the results correctly', {
  called_count <<- 0
  input_matrix <<- matrix(NA, nrow = 5, ncol = 3,
                          dimnames = list(NULL, c('rumination', 'happiness', 'activity')))
  input_matrix[, 'rumination'] <<- c(1, 0, 0, 0, 0)
  input_matrix[, 'happiness'] <<- c(0, 1, 0, 0, 0)
  input_matrix[, 'activity'] <<- c(0, 0, 1, 0, 0)
  with_mock(
    `autovarCore:::column_is_invalid` = function(...) {
      called_count <<- called_count + 1
      expect_equal(list(...), list(input_matrix[, called_count]))
      TRUE
    },
    expect_equal(autovarCore:::invalid_mask(input_matrix), 7)
  )
  expect_equal(called_count, 3)
  rm(list = c('called_count',
              'input_matrix'), pos = '.GlobalEnv')
})

test_that('invalid_mask calls column_is_invalid correctly', {
  called_count <<- 0
  one_column <<- c(0, 0, 0, 0, 1, 0, 0, 0)
  input_matrix <- matrix(NA, nrow = 8, ncol = 1, dimnames = list(NULL, 'activity'))
  input_matrix[, 'activity'] <- one_column
  with_mock(
    `autovarCore:::column_is_invalid` = function(...) {
      called_count <<- called_count + 1
      expect_equal(list(...), list(one_column))
      TRUE
    },
    expect_equal(autovarCore:::invalid_mask(input_matrix), 1)
  )
  expect_equal(called_count, 1)
  rm(list = c('called_count',
              'one_column'), pos = '.GlobalEnv')
})


test_that('column_is_invalid returns the correct result', {
  outlier_column <- c(0, 0, 0, 0, 1, 0)
  expect_false(autovarCore:::column_is_invalid(outlier_column))
  outlier_column <- rep.int(0, times = 90)
  expect_true(autovarCore:::column_is_invalid(outlier_column))
})
