context('apply_ln_transformation')

test_accuracy <- 0.0001


test_that('apply_ln_transformation calls ln_column for each column', {
  ln_column_count <<- 0
  column_names <- c('rumination', 'happiness', 'activity')
  input_matrix <- matrix(1,
                         ncol = 3,
                         nrow = 40,
                         dimnames = list(NULL, column_names))
  expected_result <- matrix(2,
                            ncol = 3,
                            nrow = 40,
                            dimnames = list(NULL, column_names))
  with_mock(
    `autovarCore:::ln_column` = function(...) {
      ln_column_count <<- ln_column_count + 1
      2
    },
    expect_equal(autovarCore:::apply_ln_transformation(input_matrix),
                 expected_result)
  )
  expect_equal(ln_column_count, 3)
  # test again without having named columns
  ln_column_count <<- 0
  input_matrix <- matrix(1,
                         ncol = 3,
                         nrow = 40)
  expected_result <- matrix(2,
                            ncol = 3,
                            nrow = 40)
  with_mock(
    `autovarCore:::ln_column` = function(...) {
      ln_column_count <<- ln_column_count + 1
      2
    },
    expect_equal(autovarCore:::apply_ln_transformation(input_matrix),
                 expected_result)
  )
  expect_equal(ln_column_count, 3)
  rm(list = 'ln_column_count', pos = '.GlobalEnv')
})

test_that('apply_ln_transformation works with different offets for different columns', {
  column_names <- c('rumination', 'happiness')
  input_matrix <- matrix(NA, ncol = 2, nrow = 5, dimnames = list(NULL, column_names))
  input_matrix[, 1] <- c(7, 0.5, 10, 349.23, 33.42)
  input_matrix[, 2] <- c(7, -23.47, 10, 349.23, 33.42)
  expected_result <- matrix(NA, ncol = 2, nrow = 5, dimnames = list(NULL, column_names))
  expected_result[, 1] <- log(input_matrix[, 1] + 0.5)
  expected_result[, 2] <- log(input_matrix[, 2] + 24.47)
  expect_less_than(sum(abs(autovarCore:::apply_ln_transformation(input_matrix) -
                           expected_result)), test_accuracy)
})


test_that('ln_column returns the actual log', {
  expect_less_than(abs(autovarCore:::ln_column(1)), test_accuracy)
  expect_less_than(abs(autovarCore:::ln_column(exp(1)) - 1), test_accuracy)
  input_data <- c(7, 2, 10, 349.23, 33.42)
  expected_result <- log(input_data)
  expect_less_than(sum(abs(autovarCore:::ln_column(input_data) - expected_result)),
                   test_accuracy)
})

test_that('ln_column works with some NA values', {
  input_data <- c(NA, 2, NA, 349.23, 33.42)
  expected_result <- log(input_data)
  expect_less_than(sum(abs(autovarCore:::ln_column(input_data) - expected_result), na.rm = TRUE),
                   test_accuracy)
})

test_that('ln_column works for values less than 1', {
  expect_less_than(abs(autovarCore:::ln_column(0)), test_accuracy)
  input_data <- c(7, 0.5, 10, 349.23, 33.42)
  expected_result <- log(input_data + 0.5)
  expect_less_than(sum(abs(autovarCore:::ln_column(input_data) - expected_result)),
                   test_accuracy)
})

test_that('ln_column works for negative numbers', {
  expect_less_than(abs(autovarCore:::ln_column(-342)), test_accuracy)
  input_data <- c(7, -23.47, 10, 349.23, 33.42)
  expected_result <- log(input_data + 24.47)
  expect_less_than(sum(abs(autovarCore:::ln_column(input_data) - expected_result)),
                   test_accuracy)
})
