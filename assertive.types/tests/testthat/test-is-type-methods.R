test_that("test.is_class.lm_numeric_raster.returns_true", {
  x <- c("lm", "numeric", "raster")
  expected <- c(TRUE, TRUE, FALSE)
  names(expected) <- x
  expect_equal(is_class(x), expected)
})
