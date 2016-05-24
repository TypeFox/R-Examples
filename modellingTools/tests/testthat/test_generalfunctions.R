require(modellingTools, quietly = TRUE, warn.conflicts = FALSE)
require(dplyr, quietly = TRUE, warn.conflicts = FALSE)
context("General Functions")

x <- data_frame(v1 = c(1,2,3),
                v2 = factor(v1),
                v3 = c("1","2","3"))

test_that("column_vector gives appropriate content", {
  expect_equal(column_vector(x,"v1"),c(1,2,3))
  expect_equal(column_vector(x,"v2"),factor(c(1,2,3)))
  expect_equal(column_vector(x,"v3"),c("1","2","3"))
})

test_that("column_vector gives appropriate classes", {
  expect_equal(class(column_vector(x,"v1")),"numeric")
  expect_equal(class(column_vector(x,"v2")),"factor")
  expect_equal(class(column_vector(x,"v3")),"character")
})


