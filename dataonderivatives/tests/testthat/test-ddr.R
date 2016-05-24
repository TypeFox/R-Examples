context("DDR")

test_that('DDR URL scheme still valid', {
  skip_on_cran()
  for (asset_class in names(DDR_ASSET_CLASSES)) {
    expect_true(httr::url_ok(ddr_url(lubridate::ymd(20150430), asset_class)))
    expect_false(httr::url_ok(ddr_url(lubridate::today() + lubridate::days(2),
      asset_class)))
  }
})

test_that("DDR zip can be downloaded", {
  skip_on_cran()
  expect_equal(download_ddr_zip(lubridate::ymd(20150430), "IR"), 0)
  expect_equal(download_ddr_zip(lubridate::ymd(20050430), "IR"), -1)
})

test_that("DDR file parsed correctly", {
  skip_on_cran()
  for (asset_class in names(DDR_ASSET_CLASSES)) {
    res <- get_ddr_data(lubridate::ymd(20150430), asset_class, TRUE)
    expect_equal(nrow(readr::problems(res)), 0)
  }
})
