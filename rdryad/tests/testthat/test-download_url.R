context("download_url")

test_that("download_url basic functionality works", {
  skip_on_cran()

  aa <- download_url(id = '10255/dryad.1759')
  bb <- download_url(id = '10255/dryad.102551')

  expect_is(aa, "character")
  expect_is(bb, "character")
  expect_true(grepl("datadryad.org/bitstream", aa))
  expect_true(grepl("datadryad.org/bitstream", bb))
  expect_true(grepl("dryad\\.1759", aa))
  expect_true(grepl("dryad\\.102551", bb))
})


test_that("download_url fails well", {
  skip_on_cran()

  expect_error(download_url(id = "10255/dryad.1664"), "No output from search")
})
