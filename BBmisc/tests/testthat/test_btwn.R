context("btwn")

test_that("btwn", {
  y = c(-1L,5L,Inf)
  expect_equal(1L:3L %btwn% y, c(TRUE, TRUE, TRUE))
  expect_equal(-2L:-1L %btwn% y, c(FALSE,TRUE))
  y = 5L
  expect_equal(5L %btwn% y, TRUE)
  expect_equal(1L %btwn% y, FALSE)
})
