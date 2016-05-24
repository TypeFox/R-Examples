context('config')

# Configuration and defaults

test_that('default_autovar_params returns a list of expected commands', {
  default_params <- autovarCore:::default_autovar_params()
  expect_equal(class(default_params), 'list')
  expect_more_than(length(names(default_params)), 4)
})

test_that('supported_test_names returns a character vector', {
  supported_test_names <- autovarCore:::supported_test_names()
  expect_equal(class(supported_test_names), 'character')
  expect_more_than(length(supported_test_names), 3)
})

test_that('run_test accepts exactly the correct test_names', {
  test_names <- autovarCore:::supported_test_names()
  for (test_name in test_names)
    expect_false(is.null(autovarCore:::run_test(test_name)))
  expect_error(autovarCore:::run_test(), 'argument .* is missing')
  expect_error(autovarCore:::run_test('asdf'), 'Unknown test')
})

test_that('supported_criteria returns a character vector', {
  supported_criteria <- autovarCore:::supported_criteria()
  expect_equal(class(supported_criteria), 'character')
  expect_more_than(length(supported_criteria), 1)
})

test_that('p_level_for_trend_significance returns a single float', {
  p_level_for_trend_significance <- autovarCore:::p_level_for_trend_significance()
  expect_equal(class(p_level_for_trend_significance), 'numeric')
  expect_equal(length(p_level_for_trend_significance), 1)
})

test_that('std_factor_for_normal_outliers returns a single float', {
  std_factor_for_normal_outliers <- autovarCore:::std_factor_for_normal_outliers()
  expect_equal(class(std_factor_for_normal_outliers), 'numeric')
  expect_equal(length(std_factor_for_normal_outliers), 1)
})

test_that('std_factor_for_squared_outliers returns a single float', {
  std_factor_for_squared_outliers <- autovarCore:::std_factor_for_squared_outliers()
  expect_equal(class(std_factor_for_squared_outliers), 'numeric')
  expect_equal(length(std_factor_for_squared_outliers), 1)
})
