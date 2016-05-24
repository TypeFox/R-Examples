context('impute_datamatrix')

test_imputation_iterations <- 3

testdata_matrix_with_missings <- function() {
  data_matrix <- matrix(nrow = 40, ncol = 3)
  data_matrix[, ] <- runif(ncol(data_matrix) * nrow(data_matrix), 1, nrow(data_matrix))
  while (sum(is.na(data_matrix)) == 0)
    data_matrix[as.logical(round(runif(ncol(data_matrix) * nrow(data_matrix), -0.3, 0.7)))] <- NA
  colnames(data_matrix) <- c('rumination', 'happiness', 'activity')
  data_matrix
}

testdata_matrix_without_missings <- function() {
  data_matrix <- matrix(nrow = 40, ncol = 3)
  data_matrix[, ] <- runif(ncol(data_matrix) * nrow(data_matrix), 1, nrow(data_matrix))
  colnames(data_matrix) <- c('rumination', 'happiness', 'activity')
  data_matrix
}


test_that('impute_datamatrix returns the original matrix if there are no missings', {
  input_matrix <- testdata_matrix_without_missings()
  expect_equal(autovarCore:::impute_datamatrix(input_matrix, 1, test_imputation_iterations),
               input_matrix)
})

test_that('impute_datamatrix returns a matrix of the correct dimensions', {
  input_matrix <- testdata_matrix_with_missings()
  imputed_matrix <- autovarCore:::impute_datamatrix(input_matrix, 1, test_imputation_iterations)
  expect_equal(nrow(imputed_matrix), nrow(input_matrix))
  expect_equal(ncol(imputed_matrix), ncol(input_matrix))
})

test_that('impute_datamatrix works for a column with just one value not NA', {
  input_matrix <- testdata_matrix_with_missings()
  input_matrix[2, 3] <- 27
  # Assert that there actually are missings:
  expect_more_than(sum(is.na(input_matrix)), 0)
  expect_equal(sum(is.na(autovarCore:::impute_datamatrix(input_matrix,
                                                         3,
                                                         test_imputation_iterations))),
               0)
})

test_that('impute_datamatrix imputes all missing values', {
  input_matrix <- testdata_matrix_with_missings()
  # Assert that there actually are missings:
  expect_more_than(sum(is.na(input_matrix)), 0)
  expect_equal(sum(is.na(autovarCore:::impute_datamatrix(input_matrix,
                                                         1,
                                                         1))),
               0)
  expect_equal(sum(is.na(autovarCore:::impute_datamatrix(input_matrix,
                                                         3,
                                                         test_imputation_iterations))),
               0)
})


test_that('impute_datamatrix does not change non missings', {
  input_matrix <- testdata_matrix_with_missings()
  expect_less_than(sum(abs(input_matrix -
                           autovarCore:::impute_datamatrix(input_matrix,
                                                           1,
                                                           test_imputation_iterations)),
                       na.rm = TRUE), 0.0001)
})

test_that('impute_datamatrix is able to handle constant columns', {
  input_matrix <- testdata_matrix_with_missings()
  input_matrix[, 'happiness'] <- 20
  imputed_matrix <- autovarCore:::impute_datamatrix(input_matrix,
                                                    1,
                                                    test_imputation_iterations)
  expect_less_than(sum(abs(input_matrix - imputed_matrix), na.rm = TRUE), 0.0001)
  expect_equal(sum(is.na(imputed_matrix)), 0)
  input_matrix[c(1, 9, 23, 36), 'happiness'] <- NA
  imputed_matrix <- autovarCore:::impute_datamatrix(input_matrix,
                                                    1,
                                                    test_imputation_iterations)
  expect_less_than(sum(abs(input_matrix - imputed_matrix), na.rm = TRUE), 0.0001)
  expect_equal(sum(is.na(imputed_matrix)), 0)
})

test_that('impute_datamatrix is able to handle some columns not having missing values', {
  input_matrix <- testdata_matrix_with_missings()
  without_missing_data <- testdata_matrix_without_missings()
  input_matrix[, 'rumination'] <- without_missing_data[, 'rumination']
  imputed_matrix <- autovarCore:::impute_datamatrix(input_matrix,
                                                    1,
                                                    test_imputation_iterations)
  expect_less_than(sum(abs(input_matrix - imputed_matrix), na.rm = TRUE), 0.0001)
  expect_equal(sum(is.na(imputed_matrix)), 0)
  input_matrix[, 'activity'] <- without_missing_data[, 'activity']
  imputed_matrix <- autovarCore:::impute_datamatrix(input_matrix,
                                                    1,
                                                    test_imputation_iterations)
  expect_less_than(sum(abs(input_matrix - imputed_matrix), na.rm = TRUE), 0.0001)
  expect_equal(sum(is.na(imputed_matrix)), 0)
})

test_that('impute_datamatrix calls amelia the specified amount of times', {
  expected_result <- matrix(1,
                            nrow = 40,
                            ncol = 3,
                            dimnames = list(NULL, c('rumination', 'happiness', 'activity')))
  amelia_counter <<- 0
  with_mock(
    `Amelia::amelia` = function(...) {
      amelia_counter <<- amelia_counter + 1
      list(imputations = 5)
    },
    expect_equal(autovarCore:::impute_datamatrix(testdata_matrix_with_missings(),
                                                 1, 5),
                                                 expected_result)
  )
  expect_equal(amelia_counter, 5)
  amelia_counter <<- 0
  with_mock(
    `Amelia::amelia` = function(...) {
      amelia_counter <<- amelia_counter + 1
      list(imputations = 5)
    },
    expect_equal(autovarCore:::impute_datamatrix(testdata_matrix_with_missings(),
                                                 1, 1),
                 expected_result)
  )
  expect_equal(amelia_counter, 1)
  rm(list = 'amelia_counter', pos = '.GlobalEnv')
})

test_that('has_missings detects missings correctly', {
  expect_equal(autovarCore:::has_missings(testdata_matrix_with_missings()),
               TRUE)
  expect_equal(autovarCore:::has_missings(testdata_matrix_without_missings()),
               FALSE)
})
