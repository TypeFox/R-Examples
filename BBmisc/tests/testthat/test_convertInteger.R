context("convert ints")

test_that("convertInteger", {
  expect_true(identical(convertInteger(1), 1L))
  expect_true(identical(convertInteger(1L), 1L))
  expect_true(identical(convertInteger(c(1,4)), c(1, 4)))
  expect_true(identical(convertInteger("a"), "a"))
  expect_true(identical(convertInteger(NA), as.integer(NA)))
  expect_true(identical(convertInteger(as.integer(NA)), as.integer(NA)))
  expect_true(identical(convertInteger(as.numeric(NA)), as.integer(NA)))
  expect_true(identical(convertInteger(c(1, NA)), c(1, NA)))
})


test_that("convertIntegers", {
  expect_true(identical(convertIntegers(1), 1L))
  expect_true(identical(convertIntegers(1L), 1L))
  expect_true(identical(convertIntegers(c(1,4)), c(1L, 4L)))
  expect_true(identical(convertIntegers("a"), "a"))
  expect_true(identical(convertIntegers(NA), as.integer(NA)))
  expect_true(identical(convertIntegers(c(NA, NA)), as.integer(c(NA, NA))))
  expect_true(identical(convertIntegers(as.integer(c(NA, NA))), as.integer(c(NA, NA))))
  expect_true(identical(convertIntegers(as.numeric(c(NA, NA))), as.integer(c(NA, NA))))
  expect_true(identical(convertIntegers(c(1, NA)), as.integer(c(1, NA))))  
  expect_true(identical(convertIntegers(c()), integer()))
  expect_true(identical(convertIntegers(c(x = 1, y = 4)), c(x = 1L, y = 4L)))
})
