### Assignment : task2 ###

context("task2")

test_that("Marking task2", {
  expect_true(exists("task2"), info = "task2() does not exist.")
  expect_is(task2, "function", info = "task2 is not a function.")
  expect_function_arguments(task2, c("vector"))
})
