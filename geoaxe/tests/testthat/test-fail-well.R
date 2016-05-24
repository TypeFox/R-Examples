context("geoaxe fails well")

test_that("fails wells", {
  expect_error(chop(1), "method not implemented for numeric")
  expect_error(chop(mtcars), "method not implemented for data.frame")
  expect_error(chop(matrix()), "method not implemented for matrix")
  expect_error(chop(NULL), "method not implemented for NULL")
  expect_error(chop("asdfafsf"), "input must be WKT or GeoJSON")
})
