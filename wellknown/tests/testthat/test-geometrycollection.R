context("geometrycollection")

test_that("geometrycollection works", {
  ## empty gc
  empty <- geometrycollection("empty")
  expect_is(empty, "character")
  expect_equal(empty, "GEOMETRYCOLLECTION EMPTY")

  ## returns self
  self <- geometrycollection("GEOMETRYCOLLECTION(POINT(4 6), LINESTRING(4 6, 7 10))")
  expect_is(self, "character")
  expect_equal(self, "GEOMETRYCOLLECTION(POINT(4 6), LINESTRING(4 6, 7 10))")

  ## from a point
  pt <- geometrycollection(point(-116.4, 45.2))
  expect_is(pt, "character")
  expect_equal(pt, "GEOMETRYCOLLECTION (POINT (-116.4000000000000057 45.2000000000000028))")

  ## from many inputs
  many <- geometrycollection(point(-116.4, 45.2, fmt = 2),
    linestring("LINESTRING (-116.4 45.2, -118.0 47.0)"),
    circularstring(list(c(1, 5), c(6, 2), c(7, 3)), fmt = 2)
  )
  many_match <- "GEOMETRYCOLLECTION (POINT (-116.40 45.20), LINESTRING (-116.4 45.2, -118.0 47.0), CIRCULARSTRING (1.00 5.00, 6.00 2.00, 7.00 3.00))"
  expect_is(many, "character")
  expect_equal(many, many_match)
})

test_that("geometrycollection fails correctly", {
  expect_error(geometrycollection(-116.4), "no applicable method")
  expect_error(geometrycollection(), "no applicable method")
  expect_error(geometrycollection(mtcars), "no applicable method")
  expect_error(geometrycollection("POINT(5)"), "All inputs must be WKT strings")
})
