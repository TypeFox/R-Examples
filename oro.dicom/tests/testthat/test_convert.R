context("Numeric conversions")

test_that("Convert integer to different base", {
  expect_equal(dec2base(23, 2), "10111")
  expect_equal(dec2base(2, 2), "10")
})

test_that("Convert integer to hexadecimal base", {
  expect_equal(dec2hex(1), "1")
  expect_equal(dec2hex(2), "2")
  expect_equal(dec2hex(3), "3")
  expect_equal(dec2hex(4), "4")
  expect_equal(dec2hex(5), "5")
  expect_equal(dec2hex(6), "6")
  expect_equal(dec2hex(7), "7")
  expect_equal(dec2hex(8), "8")
  expect_equal(dec2hex(9), "9")
  expect_equal(dec2hex(10), "A")
  expect_equal(dec2hex(11), "B")
  expect_equal(dec2hex(12), "C")
  expect_equal(dec2hex(13), "D")
  expect_equal(dec2hex(14), "E")
})

test_that("Convert strings to dates/times", {
  expect_is(str2date("19930822"), "character")
  expect_equal(str2date("19930822"), "22 Aug 1993")
  expect_is(str2time("112308")$txt, "character")
  expect_equal(str2time("112308")$txt, "11:23:08.00000")
  expect_is(str2time("112308")$time, "numeric")
  expect_equal(str2time("112308")$time, 40988)
})
  
