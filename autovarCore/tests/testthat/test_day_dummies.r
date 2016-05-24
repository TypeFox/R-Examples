context('day_dummies')

test_that('day_dummies returns NULL for 0 measurements per day or for one row', {
  expect_null(autovarCore:::day_dummies(40, 0))
  expect_null(autovarCore:::day_dummies(1, 3))
  expect_null(autovarCore:::day_dummies(1, 1))
})

test_that('day_dummies calls seasonal_dummy_columns correctly', {
  called_count <<- 0
  calling_args <<- NULL
  with_mock(
    `autovarCore:::seasonal_dummy_columns` = function(...) {
      called_count <<- called_count + 1
      calling_args <<- list(...)
      'testresult'
    },
    expect_equal(autovarCore:::day_dummies(15, 5),
                 'testresult')
  )
  expect_equal(called_count, 1)
  expect_equal(calling_args, list(out_length = 15, period = 7, repetitions = 5, dummy_name_prefix = 'day_'))
  rm(list = c('called_count', 'calling_args'), pos = '.GlobalEnv')
})

test_that('day_dummies returns the correct result', {
  column_names <- c('day_1', 'day_2')
  expected_result <- matrix(NA,
                            nrow = 9,
                            ncol = 2,
                            dimnames = list(NULL, column_names))
  expected_result[, 1] <- c(1, 1, 1, 1, 0, 0, 0, 0, 0)
  expected_result[, 2] <- c(0, 0, 0, 0, 1, 1, 1, 1, 0)
  expect_equal(autovarCore:::day_dummies(9, 4),
               expected_result)
  column_names <- c('day_1', 'day_2', 'day_3', 'day_4', 'day_5', 'day_6')
  expected_result <- matrix(NA,
                            nrow = 9,
                            ncol = 6,
                            dimnames = list(NULL, column_names))
  expected_result[, 1] <- c(1, 0, 0, 0, 0, 0, 0, 1, 0)
  expected_result[, 2] <- c(0, 1, 0, 0, 0, 0, 0, 0, 1)
  expected_result[, 3] <- c(0, 0, 1, 0, 0, 0, 0, 0, 0)
  expected_result[, 4] <- c(0, 0, 0, 1, 0, 0, 0, 0, 0)
  expected_result[, 5] <- c(0, 0, 0, 0, 1, 0, 0, 0, 0)
  expected_result[, 6] <- c(0, 0, 0, 0, 0, 1, 0, 0, 0)
  expect_equal(autovarCore:::day_dummies(9, 1),
               expected_result)
})
