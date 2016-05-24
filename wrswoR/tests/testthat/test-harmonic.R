context("harmonic")

test_that("values around length of lookup table", {
  i <- length(.harmonic.series)
  x <- sapply(seq.int(i - 5, i + 5), .harmonic)
  expect_true(all(diff(x) > 0))
})
