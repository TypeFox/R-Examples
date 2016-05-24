context("lookahead")

test_that(
  "lookahead wraps in lookahead token",
  {
    expected <- as.regex("(?=foo)")
    actual <- lookahead("foo")
    expect_equal(actual, expected)
  }
)

context("lookbehind")

test_that(
  "lookbehind wraps in lookbehind token",
  {
    expected <- as.regex("(?<=foo)")
    actual <- lookbehind("foo")
    expect_equal(actual, expected)
  }
)

context("negative_lookahead")

test_that(
  "negative_lookahead wraps in negative lookahead token",
  {
    expected <- as.regex("(?!foo)")
    actual <- negative_lookahead("foo")
    expect_equal(actual, expected)
  }
)

context("negative_lookbehind")

test_that(
  "negative_lookbehind wraps in negative lookbehind token",
  {
    expected <- as.regex("(?<!foo)")
    actual <- negative_lookbehind("foo")
    expect_equal(actual, expected)
  }
)
