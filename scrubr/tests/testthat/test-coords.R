context("coord_* functions")

df <- sample_data_1

test_that("coord_* passing lat/long vars works", {
  skip_on_cran()

  # 1 name input
  lon_name <- "decimalLongitude"
  names(df)[2] <- lon_name
  bb <- suppressMessages(dframe(df) %>% coord_incomplete(lon = lon_name))

  expect_is(bb, "data.frame")
  expect_is(bb, "dframe")
  expect_equal(names(df)[2], lon_name)
  expect_equal(names(bb)[2], "longitude")

  # both names input
  lon_name <- "x"
  lat_name <- "y"
  names(df)[2] <- lon_name
  names(df)[3] <- lat_name
  bb <- suppressMessages(dframe(df) %>% coord_incomplete(lat_name, lon_name))

  expect_is(bb, "data.frame")
  expect_is(bb, "dframe")
  expect_equal(names(df)[2], lon_name)
  expect_equal(names(df)[3], lat_name)
  expect_equal(names(bb)[2], "x")
  expect_equal(names(bb)[3], "y")
})
