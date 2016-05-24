context("char_class")

test_that(
  "char_class wraps in class token",
  {
    expected <- as.regex("[abc]")
    actual <- char_class("a", "b", "c")
    expect_equal(actual, expected)
  }
)

context("negated_char_class")

test_that(
  "negated_char_class wraps in class token with ^",
  {
    expected <- as.regex("[^abc]")
    actual <- negated_char_class("a", "b", "c")
    expect_equal(actual, expected)
  }
)
