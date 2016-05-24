context("properties")

str <- "POINT (-116.4000000000000057 45.2000000000000028)"
x <- wkt2geojson(str)
a <- properties(x, style=list(color = "red"))

test_that("set propertiers works", {
  expect_is(a, "geojson")
  expect_match(a$properties$style$color, "red")
  expect_null(a$properties$popup)
})

test_that("properties deals with bad input well", {
  expect_error(properties(x, style=list()), "needs a non-empty list")
  expect_error(properties(x, style=NULL), "supply a list of named options")
  expect_error(properties(x, popup=NULL), "supply a list of named options")
})
