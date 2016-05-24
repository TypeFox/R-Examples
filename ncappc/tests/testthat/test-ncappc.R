context("Check ncappc")

test_that("Check if ncappc is available",{
  expect_true(exists("ncappc"), TRUE)
  source("ncappcexamples/example1.R")
  expect_equal(sum(TIME),sum(TIME.positive))
  expect_equal(sum(CONC),sum(CONC.positive))
})

