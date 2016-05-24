context("BSEF")
test_that('BSEF API accesible', {
  skip_on_cran()
  empty_data <- download_bsef_data_single(lubridate::ymd(20140918), 'IR')
  expect_false(identical(empty_data, dplyr::data_frame()))
  empty_data <- get_bsef_data(lubridate::ymd(20140920))
  expect_true(identical(empty_data, dplyr::data_frame()))
})
