test_that("get_data returns a data frame", {

  key <- "ICP.M.DE.N.000000+XEF000.4.ANR"
  filter <- list(lastNObservations = 12)

  hicp <- get_data(key, filter)

  expect_equal(class(hicp), c("tbl_df", "tbl", "data.frame"))
})

test_that("get_dataflows returns a data frame", {
  expect_equal(class(get_dataflows()), c("tbl_df", "tbl", "data.frame"))
})

test_that("get_dimensions returns a list of data frames", {
  dims <- get_dimensions("ICP.M.DE.N.000000+XEF000.4.ANR")
  expect_equal(class(dims), "list")

  dim_classes <- vapply(dims, function(x) class(x) == "data.frame", logical(1))
  expect_true(all(dim_classes))
})

test_that("malformed series key returns 404 error", {
  key <- "ICP.M.DE.N.000000+XEF000.4.ANRs"
  expect_error(get_data(key), regexp = "404")
  expect_error(get_dimensions(key), regexp = "404")
  expect_error(get_description(key), regexp = "404")
})

test_that("get_description returns a character vector", {
  key <- "STS.A.AT+DE.N.UNEH.RTT000.4.AV3"
  desc <- get_description(key)
  expect_true(is.character(desc))

  desc_lengths <- vapply(desc, nchar, numeric(1))
  expect_true(all(desc > 0))
})
