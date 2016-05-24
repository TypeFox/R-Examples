context("datetime")

test_that(
  "datetime returns a correct regex for '%Y-%m-%d %H:%M:%S'",
  {
    expected <- as.regex("(?:[0-9]{4}-(?:0[1-9]|1[0-2])-(?:0[1-9]|[12][0-9]|3[01]) (?:[01][0-9]|2[0-3]):[0-5][0-9]:(?:[0-5][0-9]|6[01]))")
    actual <- datetime("%Y-%m-%d %H:%M:%S")
    expect_equal(actual, expected)
  }
)

