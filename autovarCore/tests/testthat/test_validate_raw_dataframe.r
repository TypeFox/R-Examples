context('validate_raw_dataframe')

testdata_raw_dataframe <- function() {
  data.frame(id = rep(1, times = 5),
             tijdstip = c("1", "3", "5", "6", "7"),
             home = c(1, 0, 1, NA, 1),
             stringsAsFactors = FALSE)
}

testdata_data_matrix <- function() {
  data_matrix <- matrix(ncol = 3, nrow = 5)
  data_matrix[, 1] <- 1
  data_matrix[, 2] <- c(1, 3, 5, 6, 7)
  data_matrix[, 3] <- c(1, 0, 1, NA, 1)
  colnames(data_matrix) <- c('id', 'tijdstip', 'home')
  data_matrix
}


test_that('validate_raw_dataframe does not accept a null data frame', {
  expect_error(autovarCore:::validate_raw_dataframe(),
               "is missing, with no default")
  expect_error(autovarCore:::validate_raw_dataframe(NULL),
               "Given param cannot be NULL")
})

test_that('validate_raw_dataframe requires class data.frame', {
  expect_error(autovarCore:::validate_raw_dataframe(list(a = 2)),
               "Param class should be: data.frame")
})

test_that('validate_raw_dataframe requires at least one row', {
  expect_error(autovarCore:::validate_raw_dataframe(data.frame()),
               "The number of rows in the data frame is below the minimum of 1")
})

test_that('validate_raw_dataframe returns a numeric matrix with labels', {
  expected_result <- testdata_data_matrix()
  expect_equal(autovarCore:::validate_raw_dataframe(testdata_raw_dataframe()),
               expected_result)
})
