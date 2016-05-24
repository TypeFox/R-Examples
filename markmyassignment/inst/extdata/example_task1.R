### Assignment : task1 ###

context("task1")

test_that("Marking of task 1", {
  expect_true(exists("task1"), info = "task1() does not exist.")
  expect_is(task1, "numeric", info = "task1 is not return a numeric value")
  expect_equal(task1, c(pi, exp(1)), info = "task1 contains erroneous values")
})
