context("symdiff")

test_that("symdiff", {
  expect_equal(symdiff(c(1, 2), 1), 2)
  expect_equal(symdiff(c(1, 2), numeric(0)), c(1, 2))
  expect_equal(symdiff("a", "b"), c("a", "b"))
  expect_equal(symdiff("a", "a"), character(0))
})


