context("Variability Helper Functions")

test_that("rmssd estimates correctly", {
  expect_that(rmssd(1:4), equals(1))
})

test_that("rmssd by ID estimates correctly", {
  expect_that(rmssd_id(c(1:3, 1, 3, 5), rep(1:2, each = 3), FALSE), is_equivalent_to(array(c(1, 2))))
})

test_that("sd by ID estimates correctly", {
  expect_that(sd_id(c(1:3, 1, 3, 5), rep(1:2, each = 3), FALSE), is_equivalent_to(array(c(1, 2))))
})

test_that("rolling_diff estimates correctly", {
  expect_that(rolling_diff(1:4, 4), is_equivalent_to(3))
})

test_that("rolling_diff by ID estimates correctly", {
  expect_that(rolling_diff_id(c(1:4, 1, 3, 5, 7), rep(1:2, each = 4), FALSE, 4), is_equivalent_to(array(c(3, 6))))
})

test_that("empirical_pvalue returns correct p-values", {
  expect_that(empirical_pvalue(c(-1, 1:3))["p-value"], is_equivalent_to(.5))
})
