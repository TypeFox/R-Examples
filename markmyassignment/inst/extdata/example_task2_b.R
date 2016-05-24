### Assignment : task2 ###

context("task2")

test_that("Mark even more on task2", {
  expect_is(task2(5:10), "integer", info = "task2 don't return an integer.")
  expect_equal(length(task2(5:10)), 1, info = "task2 do not return one value.")
  expect_equal(task2(vector=5:10), 15, info = "task2(vector=5:10) don't return 15")
  expect_equal(task2(vector=1:5), 6, info = "task2(vector=1:5) don't return 6")
})
