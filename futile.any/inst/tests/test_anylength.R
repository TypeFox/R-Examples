context("anylength")

test_that("vector has correct length", {
  a <- c(1,2,3)
  expect_that(anylength(a), equals(3))
})

