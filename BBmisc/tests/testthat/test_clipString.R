context("clipString")

test_that("clipString", {
  expect_equal(clipString("abcdef", 7), "abcdef")
  expect_equal(clipString("abcdef", 6), "abcdef")
  expect_equal(clipString("abcdef", 5), "ab...")
  expect_error(clipString("abcdef", 2))

  expect_equal(clipString(NA_character_, 5), NA_character_)
  expect_equal(clipString(c("aaaaa", NA, "aa"), 4), c("a...", NA, "aa"))
})
