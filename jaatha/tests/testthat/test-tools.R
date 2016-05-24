context("Tools")

test_that("it checks for packages", {
  skip_on_cran()
  expect_true(require_package("jaatha"))
  expect_error(require_package("2l3ihjrpaiwhf"))
})


test_that("it samples seeds", {
  set.seed(17)
  seeds <- sample_seed(10)
  expect_equal(length(seeds), 10)
  expect_equal(length(unique(seeds)), 10)
  
  set.seed(17)
  seeds2 <- sample_seed(10)
  expect_equal(seeds2, seeds)
  
  seeds <- sample_seed(7)
  expect_equal(length(seeds), 7)
})
