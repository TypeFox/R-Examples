context("metric info")

test_that("a string has positive dimensions", {
  value <- str_extents("Hello World!")
  expect_true(all( value > 0 ))
})

