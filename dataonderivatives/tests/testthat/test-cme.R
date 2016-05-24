context("CME SDR")

test_that("CME SDR URL works",  {
  skip_on_cran()
  expect_equal(cme_ftp_url(lubridate::ymd(20150301), "FX"),
    "ftp://ftp.cmegroup.com/sdr/fx/2015/03/RT.FX.20150301.csv.zip")
  expect_equal(cme_ftp_url(lubridate::ymd(20130702), "IR"),
    "ftp://ftp.cmegroup.com/sdr/rates/2013/07/RT.IRS.20130702.csv.zip")
  expect_equal(cme_ftp_url(lubridate::ymd(20150201), "CO"),
    "ftp://ftp.cmegroup.com/sdr/commodities/2015/02/RT.COMMODITY.20150201.csv.zip")
})

test_that("CME SDR download works",  {
  skip_on_cran()
  expect_equal(download_cme_zip(lubridate::ymd(20150301), "FX"), 0)
  expect_equal(download_cme_zip(lubridate::ymd(20050301), "FX"), -1)
})


test_that("CME SDR file parses",  {
  skip_on_cran()
  expect_true(nrow(get_cme_data(lubridate::ymd(20150301), "FX")) > 0)
  expect_true(nrow(get_cme_data(lubridate::ymd(20050301), "FX")) == 0)
})
