test_that("result of function", {
  expect_equal(class(deputies_get_data("active", 7)), "data.frame")
})

test_that("result of function", {
  expect_equal(class(deputies_get_data("inactive", 8)), "data.frame")
})

test_that("columns of table", {
  expect_equal(ncol(deputies_get_data("active", 7)), 2)
})

test_that("columns of table", {
  expect_equal(ncol(deputies_get_data("inactive", 8)), 2)
})

test_that("rows of table", {
  expect_more_than(nrow(deputies_get_data("active", 7)), 0)
})

test_that("rows of table", {
  expect_more_than(nrow(deputies_get_data("inactive", 8)), 0)
})