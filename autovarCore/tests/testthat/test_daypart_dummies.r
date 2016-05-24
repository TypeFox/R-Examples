context('daypart_dummies')

test_that('daypart_dummies returns NULL for 0 or 1 measurements per day', {
  expect_null(autovarCore:::daypart_dummies(40, 0))
  expect_null(autovarCore:::daypart_dummies(13, 1))
})

test_that('daypart_dummies calls seasonal_dummy_columns correctly', {
  called_count <<- 0
  calling_args <<- NULL
  with_mock(
    `autovarCore:::seasonal_dummy_columns` = function(...) {
      called_count <<- called_count + 1
      calling_args <<- list(...)
      'testresult'
    },
    expect_equal(autovarCore:::daypart_dummies(15, 5),
                 'testresult')
  )
  expect_equal(called_count, 1)
  expect_equal(calling_args, list(out_length = 15, period = 5, repetitions = 1, dummy_name_prefix = 'dailymeas_'))
  rm(list = c('called_count', 'calling_args'), pos = '.GlobalEnv')
})

test_that('daypart_dummies returns the correct result', {
  column_names <- c('dailymeas_1', 'dailymeas_2', 'dailymeas_3')
  expected_result <- matrix(NA,
                            nrow = 9,
                            ncol = 3,
                            dimnames = list(NULL, column_names))
  expected_result[, 1] <- c(1, 0, 0, 0, 1, 0, 0, 0, 1)
  expected_result[, 2] <- c(0, 1, 0, 0, 0, 1, 0, 0, 0)
  expected_result[, 3] <- c(0, 0, 1, 0, 0, 0, 1, 0, 0)
  expect_equal(autovarCore:::daypart_dummies(9, 4),
               expected_result)
})

test_that('seasonal_dummy_columns works with different out lengths', {
  expected_result <- matrix(NA, nrow = 4, ncol = 3, dimnames = list(NULL, c('dummy1', 'dummy2', 'dummy3')))
  expected_result[, 1] <- c(1, 0, 0, 0)
  expected_result[, 2] <- c(0, 1, 0, 0)
  expected_result[, 3] <- c(0, 0, 1, 0)
  expect_equal(autovarCore:::seasonal_dummy_columns(4, 5, 1, 'dummy'), expected_result)
  expected_result <- matrix(NA, nrow = 9, ncol = 4, dimnames = list(NULL, c('dummy1', 'dummy2', 'dummy3', 'dummy4')))
  expected_result[, 1] <- c(1, 0, 0, 0, 0, 1, 0, 0, 0)
  expected_result[, 2] <- c(0, 1, 0, 0, 0, 0, 1, 0, 0)
  expected_result[, 3] <- c(0, 0, 1, 0, 0, 0, 0, 1, 0)
  expected_result[, 4] <- c(0, 0, 0, 1, 0, 0, 0, 0, 1)
  expect_equal(autovarCore:::seasonal_dummy_columns(9, 5, 1, 'dummy'), expected_result)
})

test_that('seasonal_dummy_columns works with different periods', {
  expected_result <- matrix(NA, nrow = 6, ncol = 4, dimnames = list(NULL, c('dummy1', 'dummy2', 'dummy3', 'dummy4')))
  expected_result[, 1] <- c(1, 0, 0, 0, 0, 1)
  expected_result[, 2] <- c(0, 1, 0, 0, 0, 0)
  expected_result[, 3] <- c(0, 0, 1, 0, 0, 0)
  expected_result[, 4] <- c(0, 0, 0, 1, 0, 0)
  expect_equal(autovarCore:::seasonal_dummy_columns(6, 5, 1, 'dummy'), expected_result)
  expected_result <- matrix(NA, nrow = 6, ncol = 3, dimnames = list(NULL, c('dummy1', 'dummy2', 'dummy3')))
  expected_result[, 1] <- c(1, 0, 0, 0, 1, 0)
  expected_result[, 2] <- c(0, 1, 0, 0, 0, 1)
  expected_result[, 3] <- c(0, 0, 1, 0, 0, 0)
  expect_equal(autovarCore:::seasonal_dummy_columns(6, 4, 1, 'dummy'), expected_result)
})

test_that('seasonal_dummy_columns works with different repetitions', {
  expected_result <- matrix(NA, nrow = 7, ncol = 2, dimnames = list(NULL, c('dummy1', 'dummy2')))
  expected_result[, 1] <- c(1, 0, 0, 1, 0, 0, 1)
  expected_result[, 2] <- c(0, 1, 0, 0, 1, 0, 0)
  expect_equal(autovarCore:::seasonal_dummy_columns(7, 3, 1, 'dummy'), expected_result)
  expected_result <- matrix(NA, nrow = 7, ncol = 2, dimnames = list(NULL, c('dummy1', 'dummy2')))
  expected_result[, 1] <- c(1, 1, 0, 0, 0, 0, 1)
  expected_result[, 2] <- c(0, 0, 1, 1, 0, 0, 0)
  expect_equal(autovarCore:::seasonal_dummy_columns(7, 3, 2, 'dummy'), expected_result)
})

test_that('seasonal_dummy_columns works with different dummy name prefixes', {
  expected_result <- matrix(NA, nrow = 8, ncol = 2, dimnames = list(NULL, c('dummy1', 'dummy2')))
  expected_result[, 1] <- c(1, 1, 1, 0, 0, 0, 0, 0)
  expected_result[, 2] <- c(0, 0, 0, 1, 1, 1, 0, 0)
  expect_equal(autovarCore:::seasonal_dummy_columns(8, 3, 3, 'dummy'), expected_result)
  expected_result <- matrix(NA, nrow = 8, ncol = 2, dimnames = list(NULL, c('othername_1', 'othername_2')))
  expected_result[, 1] <- c(1, 1, 1, 0, 0, 0, 0, 0)
  expected_result[, 2] <- c(0, 0, 0, 1, 1, 1, 0, 0)
  expect_equal(autovarCore:::seasonal_dummy_columns(8, 3, 3, 'othername_'), expected_result)
})

test_that('seasonal_dummy_columns calls seasonal_dummy_column and dummy_column_names correctly', {
  seasonal_dummy_column_call_count <<- 0
  dummy_column_names_call_count <<- 0
  with_mock(
    `autovarCore:::seasonal_dummy_column` = function(...) {
      seasonal_dummy_column_call_count <<- seasonal_dummy_column_call_count + 1
      1
    },
    `autovarCore:::dummy_column_names` = function(...) {
      dummy_column_names_call_count <<- dummy_column_names_call_count + 1
      c('a', 'b', 'c')
    },
    expect_equal(autovarCore:::seasonal_dummy_columns(16, 4, 5, 'test'),
                 matrix(1, nrow = 1, ncol = 3, dimnames = list(NULL, c('a', 'b', 'c'))))
  )
  expect_equal(seasonal_dummy_column_call_count, 3)
  expect_equal(dummy_column_names_call_count, 1)
  rm(list = c('dummy_column_names_call_count', 'seasonal_dummy_column_call_count'), pos = '.GlobalEnv')
})

test_that('seasonal_dummy_column works with different lengths', {
  expect_equal(autovarCore:::seasonal_dummy_column(5, 3, 1, 0),
               c(1, 0, 0, 1, 0))
  expect_equal(autovarCore:::seasonal_dummy_column(6, 3, 1, 0),
               c(1, 0, 0, 1, 0, 0))
  expect_equal(autovarCore:::seasonal_dummy_column(7, 3, 1, 0),
               c(1, 0, 0, 1, 0, 0, 1))
  expect_equal(autovarCore:::seasonal_dummy_column(3, 3, 1, 0),
               c(1, 0, 0))
  expect_equal(autovarCore:::seasonal_dummy_column(2, 3, 1, 0),
               c(1, 0))
  expect_equal(autovarCore:::seasonal_dummy_column(1, 3, 1, 0),
               1)
})

test_that('seasonal_dummy_column works with different offsets', {
  expect_equal(autovarCore:::seasonal_dummy_column(5, 3, 1, 1),
               c(0, 1, 0, 0, 1))
  expect_equal(autovarCore:::seasonal_dummy_column(10, 5, 1, 2),
               c(0, 0, 1, 0, 0, 0, 0, 1, 0, 0))
})

test_that('seasonal_dummy_column works with different periods', {
  expect_equal(autovarCore:::seasonal_dummy_column(5, 2, 1, 0),
               c(1, 0, 1, 0, 1))
  expect_equal(autovarCore:::seasonal_dummy_column(7, 7, 1, 1),
               c(0, 1, 0, 0, 0, 0, 0))
  expect_equal(autovarCore:::seasonal_dummy_column(9, 4, 1, 0),
               c(1, 0, 0, 0, 1, 0, 0, 0, 1))
})

test_that('seasonal_dummy_column works with different repetitions', {
  expect_equal(autovarCore:::seasonal_dummy_column(9, 2, 1, 0),
               c(1, 0, 1, 0, 1, 0, 1, 0, 1))
  expect_equal(autovarCore:::seasonal_dummy_column(9, 2, 2, 1),
               c(0, 0, 1, 1, 0, 0, 1, 1, 0))
  expect_equal(autovarCore:::seasonal_dummy_column(9, 2, 3, 0),
               c(1, 1, 1, 0, 0, 0, 1, 1, 1))
})

test_that('dummy_column_names works with just one column', {
  expected_result <- 'dailymeas_1'
  expect_equal(autovarCore:::dummy_column_names(1, 'dailymeas_'), expected_result)
})

test_that('dummy_column_names works with more than just one column', {
  expected_result <- c('dailymeas1', 'dailymeas2', 'dailymeas3')
  expect_equal(autovarCore:::dummy_column_names(3, 'dailymeas'), expected_result)
})
