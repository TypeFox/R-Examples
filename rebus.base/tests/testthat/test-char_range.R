context("char_range")

test_that(
  "char_range with hi > lo creates a range, wraps in class",
  {
    expected <- as.regex("[a-c]")
    actual <- char_range("a", "c")
    expect_equal(actual, expected)
  }
)

test_that(
  "char_range with hi == lo returns lo, with a warning",
  {
    expected <- as.regex("a")
    expect_warning(
      actual <- char_range("a", "a"),
      "'lo' and 'hi' are the same value.  Return 'lo'"
    )
    expect_equal(actual, expected)
  }
)

test_that(
  "char_range with hi < lo throws an error",
  {
    expect_error(
      char_range("c", "a"),
      "'hi' is less than 'lo'"
    )
  }
)
