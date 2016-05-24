context('trend_columns')

test_that('trend_columns calls add_linear_trend and add_squared_trend', {
  add_linear_trend_counter <<- 0
  add_squared_trend_counter <<- 0
  with_mock(
    `autovarCore:::add_linear_trend` = function(...) {
      add_linear_trend_counter <<- add_linear_trend_counter + 1
      1
    },
    `autovarCore:::add_squared_trend` = function(...) {
      add_squared_trend_counter <<- add_squared_trend_counter + 1
      2
    },
    expect_equal(autovarCore:::trend_columns(17), 2)
  )
  expect_equal(add_linear_trend_counter, 1)
  expect_equal(add_squared_trend_counter, 1)
  rm(list = c('add_linear_trend_counter', 'add_squared_trend_counter'), pos = '.GlobalEnv')
})

test_that('trend_columns return a 2 column matrix with the expected number of rows', {
  expected_result <- matrix(NA, ncol = 2, nrow = 7, dimnames = list(NULL, c('index', 'index2')))
  expected_result[, 1] <- 1:7
  expected_result[, 2] <- (1:7)^2
  expect_equal(autovarCore:::trend_columns(7), expected_result)
})

test_that('add_linear_trend adds a matrix column with the correct contents', {
  expected_result <- matrix(1:10, dimnames = list(NULL, 'index'))
  expect_equal(autovarCore:::add_linear_trend(NULL, 10), expected_result)
  expected_result <- matrix(1:20, dimnames = list(NULL, 'index'))
  expect_equal(autovarCore:::add_linear_trend(NULL, 20), expected_result)
})

test_that('add_linear_trend works when the first argument is a matrix', {
  some_val_data <- c(3, 2, 5, 8, 9)
  expected_result <- matrix(NA, ncol = 2, nrow = 5, dimnames = list(NULL, c('some_val', 'index')))
  expected_result[, 'some_val'] <- some_val_data
  expected_result[, 'index'] <- 1:5
  input_matrix <- matrix(some_val_data, dimnames = list(NULL, 'some_val'))
  expect_equal(autovarCore:::add_linear_trend(input_matrix, 5), expected_result)
})

test_that('add_squared_trend adds a matrix column with the correct contents', {
  expected_result <- matrix((1:13)^2, dimnames = list(NULL, 'index2'))
  expect_equal(autovarCore:::add_squared_trend(NULL, 13), expected_result)
  expected_result <- matrix((1:40)^2, dimnames = list(NULL, 'index2'))
  expect_equal(autovarCore:::add_squared_trend(NULL, 40), expected_result)
})

test_that('add_squared_trend works when the first argument is a matrix', {
  some_val_data <- c(3, 2, 5, 8, 9)
  expected_result <- matrix(NA, ncol = 2, nrow = 5, dimnames = list(NULL, c('some_val', 'index2')))
  expected_result[, 'some_val'] <- some_val_data
  expected_result[, 'index2'] <- c(1, 4, 9, 16, 25)
  input_matrix <- matrix(some_val_data, dimnames = list(NULL, 'some_val'))
  expect_equal(autovarCore:::add_squared_trend(input_matrix, 5), expected_result)
})
