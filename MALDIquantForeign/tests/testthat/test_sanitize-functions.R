context("sanitize")


test_that("sanitize-functions", {
  s <- c("a<b", "b>c", "a&o", "abc")
  r <- c("a&lt;b", "b&gt;c", "a&amp;o", "abc")
  expect_identical(MALDIquantForeign:::.sanitize(s), r)
})
