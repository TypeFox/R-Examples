context("number_range")

test_that(
  "number_range returns a correct regex for 0 to 9999",
  {
    expected <- as.regex("(?:[0-9]{4})")
    actual <- number_range(0, 9999)
    expect_equal(actual, expected)
  }
)

test_that(
  "number_range returns a correct regex for 123 to 876",
  {
    expected <- as.regex(
      "(?:12[3-9]|1[3-9][0-9]|[2-7][0-9]{2}|8[0-6][0-9]|87[0-6])"
    )
    actual <- number_range(123, 876)
    expect_equal(actual, expected)
  }
)

test_that(
  "number_range returns a correct regex for 1010 to 9090",
  {
    expected <- as.regex(
      "(?:10[1-9][0-9]|1[1-9][0-9]{2}|[2-8][0-9]{3}|90[0-8][0-9]|9090)"
    )
    actual <- number_range(1010, 9090)
    expect_equal(actual, expected)
  }
)

test_that(
  "number_range returns a correct regex for -123 to 876",
  {
    expected <- as.regex(
      "(?:-(?:[1-9]|[1-9][0-9]|1[0-1][0-9]|12[0-3])|(?:[0-7][0-9]{2}|8[0-6][0-9]|87[0-6]))"
    )
    actual <- number_range(-123, 876)
    expect_equal(actual, expected)
  }
)
