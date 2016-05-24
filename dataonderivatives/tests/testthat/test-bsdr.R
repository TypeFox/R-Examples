context("BSDR")

test_that('BSDR API accesible', {
  skip_on_cran()
  df1 <- download_bsdr_data_single(lubridate::ymd(20150505), "CO")
  expect_is(df1, "data.frame")
  df2 <- download_bsdr_data_single(lubridate::ymd(20150504, 20150506), "CR")
  expect_is(df2, "data.frame")
  expect_true(nrow(df1) <= nrow(df2))
  df3 <- download_bsdr_data_single(lubridate::ymd(20150505), c("CO", "CR"))
  expect_is(df3, "data.frame")
  expect_true(nrow(df1) <= nrow(df3))
})
