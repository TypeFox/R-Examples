context("isValidName")

test_that("isValidName", {
  expect_true(isValidName("a"))
  expect_true(all(isValidName(c("a", "b"))))
  expect_equal(isValidName(c("a", "a")), c(TRUE, FALSE))
  expect_true(all(isValidName(c("a", "a"), unique=FALSE)))
  expect_equal(isValidName(c("x", "..1")), c(TRUE, FALSE))
})
