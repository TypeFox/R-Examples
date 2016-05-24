context("toRangeStr")

expect_range = function(x, str, ...) {
  expect_equal(toRangeStr(x, ...), str)
  expect_equal(toRangeStr(sample(x), ...), str)
  expect_equal(toRangeStr(sample(c(x, x)), ...), str)
  expect_equal(toRangeStr(sample(c(x, x)), ...), str)
}

test_that("continuous ranges", {
  x = c(1, 2, 3, 4, 5, 6)
  expect_range(x, "1 - 6")
  expect_range(x, "1_6", range.sep="_")
})

test_that("single number", {
  x = 1
  expect_range(x, "1")
})

test_that("negative numbers", {
  x = -2:4
  expect_range(x, "-2 - 4")
})

test_that("noncontinuous ranges", {
  x = c(-5, -4, -2, 0, 2, 3, 4, 7)
  expect_range(x, "-5 - -4, -2, 0, 2 - 4, 7")
})


