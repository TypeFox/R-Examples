context('assertions')

testdata_raw_dataframe <- function() {
  data.frame(id = rep(1, times = 5),
             tijdstip = c(1, 3, 5, 6, 7),
             home = c('yes', 'no', 'yes', NA, 'yes'),
             stringsAsFactors = FALSE)
}


test_that('assert_param_presence requires the presence of a param name', {
  expect_error(autovarCore:::assert_param_presence('c', c('a', 'b')),
               "c is a required parameter")
  expect_error(autovarCore:::assert_param_presence('c', 'a'),
               "c is a required parameter")
  # The statement below should not throw an error.
  expect_null(autovarCore:::assert_param_presence('a', c('a', 'b')))
  expect_null(autovarCore:::assert_param_presence('c', 'c'))
})

test_that('assert_param_class asserts that the given param is of the specified class', {
  expected_error_message <- 'Param class should be: list'
  expect_error(autovarCore:::assert_param_class(NULL, 'list'), expected_error_message)
  expect_error(autovarCore:::assert_param_class(testdata_raw_dataframe(), 'list'), expected_error_message)
  # The statement below should not throw an error.
  expect_null(autovarCore:::assert_param_class(list(), 'list'))
})

test_that('assert_param_subset throws an error if any given name is not in the allowed names', {
  expected_error_message <- "Invalid param:"
  expect_error(autovarCore:::assert_param_subset('unknown', c('a', 'b')),
               paste(expected_error_message, 'unknown'))
  expect_error(autovarCore:::assert_param_subset(c('a', 'b', 'unknown'), c('a', 'b')),
               paste(expected_error_message, 'unknown'))
  # The statements below should not throw errors.
  expect_null(autovarCore:::assert_param_subset(NULL, c('a', 'b')))
  expect_null(autovarCore:::assert_param_subset(c('a', 'b'), c('a', 'b')))
  expect_null(autovarCore:::assert_param_subset('b', c('a', 'b')))
})

test_that('assert_param_not_null throws an error when supplied a null argument', {
  expected_error_message <- 'Given param cannot be NULL'
  expect_error(autovarCore:::assert_param_not_null(NULL), expected_error_message)
  # The statements below should not throw errors.
  expect_null(autovarCore:::assert_param_not_null(1))
  expect_null(autovarCore:::assert_param_not_null('a'))
})

test_that('assert_param_integer throws an error when the given argument is not an integer', {
  expected_error_message <- 'Given param is not an integer:'
  expect_error(autovarCore:::assert_param_integer(6.7), expected_error_message)
  expect_error(autovarCore:::assert_param_integer('a'), expected_error_message)
  expect_error(autovarCore:::assert_param_integer(TRUE), expected_error_message)
  # The statements below should not throw errors.
  expect_null(autovarCore:::assert_param_integer(1))
  expect_null(autovarCore:::assert_param_integer(-7))
  expect_null(autovarCore:::assert_param_integer(0))
})

test_that('assert_param_single throws an error when the given param is not a single element', {
  expected_error_message <- 'Length of given param is not 1:'
  expect_error(autovarCore:::assert_param_single(NULL), expected_error_message)
  expect_error(autovarCore:::assert_param_single(c('a', 'b')), expected_error_message)
  # The statements below should not throw errors.
  expect_null(autovarCore:::assert_param_single(3))
  expect_null(autovarCore:::assert_param_single(-6))
  expect_null(autovarCore:::assert_param_single('b'))
  expect_null(autovarCore:::assert_param_single('ando'))
  expect_null(autovarCore:::assert_param_single(7.8))
})

test_that('assert_param_range throws an error when the given integer is not in range', {
  expected_error_message <- 'The param_name has to be an integer in range 2-10'
  expect_error(autovarCore:::assert_param_range(1, 2, 10, 'param_name'), expected_error_message)
  expect_error(autovarCore:::assert_param_range(11, 2, 10, 'param_name'), expected_error_message)
  # The statements below should not throw errors.
  expect_null(autovarCore:::assert_param_range(10, 2, 10, 'param_name'))
  expect_null(autovarCore:::assert_param_range(2, 2, 10, 'param_name'))
  expect_null(autovarCore:::assert_param_range(6, 2, 10, 'param_name'))
})

test_that('assert_param_nrow checks the number of rows in data frames', {
  expect_error(autovarCore:::assert_param_nrow(testdata_raw_dataframe(),
                                               minimum = 6),
               'The number of rows in the data frame is below the minimum of 6')
  expect_error(autovarCore:::assert_param_nrow(testdata_raw_dataframe(),
                                               maximum = 4),
               'The number of rows in the data frame is above the maximum of 4')
  expect_error(autovarCore:::assert_param_nrow(testdata_raw_dataframe(),
                                               minimum = 1,
                                               maximum = 4),
               'The number of rows in the data frame is above the maximum of 4')
  # The statements below should not throw errors.
  expect_null(autovarCore:::assert_param_nrow(testdata_raw_dataframe(),
                                              minimum = 5,
                                              maximum = 5))
  expect_null(autovarCore:::assert_param_nrow(testdata_raw_dataframe(),
                                              minimum = 3))
  expect_null(autovarCore:::assert_param_nrow(testdata_raw_dataframe(),
                                              maximum = 6))
})

