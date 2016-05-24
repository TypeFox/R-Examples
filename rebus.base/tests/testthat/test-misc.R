context("exactly")

test_that(
  "exactly wraps in start and end",
  {
    expected <- as.regex("^foo$")
    actual <- exactly("foo")
    expect_equal(actual, expected)
  }
)

context("literal")

test_that(
  "literal wraps in literal tokens",
  {
    expected <- as.regex("\\Qfoo\\E")
    actual <- literal("foo")
    expect_equal(actual, expected)
  }
)
