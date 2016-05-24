context("Check colSds and rowSds")

x <- matrix(rnorm(50), 10, 5)
dimnames(x) <- list(paste0("r", 1:10), paste0("c", 1:5))


test_that("colSds and rowSds is computed correctly", {
  expect_that(apply(x, 1, sd), equals(GMCM:::rowSds(x)))
  expect_that(apply(x, 2, sd), equals(GMCM:::colSds(x)))
})

test_that("colSds and rowSds fails as intented", {
  expect_error(GMCM:::rowSds(x[, 0]), "column")
  expect_error(GMCM:::colSds(x[0, ]), "row")
})

test_that("colSds and rowSds handles zero rows and columns correctly", {
  expect_that(GMCM:::rowSds(x[0, ]), equals(apply(x[0,], 1, sd)))
  expect_that(GMCM:::colSds(x[, 0]), equals(apply(x[,0], 2, sd)))
})


