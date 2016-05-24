context("nsr functions")

test_that("nsr works", {
  skip_on_cran()

  aa <- nsr("Pinus ponderosa", "United States")

  expect_is(aa, "data.frame")
  expect_equal(aa$species, "Pinus ponderosa")
  expect_equal(aa$native_status_sources, "usda")
})

test_that("fails well", {
  skip_on_cran()

  expect_error(nsr(), "argument \"species\" is missing")
  expect_equal(NROW(nsr(species = "adadfd", country = "United States")), 0)
})
