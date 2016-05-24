context("hasAttributes")

test_that("hasAttributes", {
  x = 1:10
  attribute.names = c("size", "importance")

  expect_false(hasAttributes(x, attribute.names))

  attr(x, "size") = length(x)
  expect_false(hasAttributes(x, attribute.names))

  attr(x, "importance") = "very important"
  expect_true(hasAttributes(x, attribute.names))
})
