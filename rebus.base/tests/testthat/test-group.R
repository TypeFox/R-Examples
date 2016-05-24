context("capture")

test_that(
  "capture wraps input with capture group token",
  {
    expected <- as.regex("(foo)")
    actual <- capture("foo")
    expect_equal(actual, expected)
  }
)

context("group")

test_that(
  "group wraps input with non-capture token",
  {
    expected <- as.regex("(?:foo)")
    actual <- group("foo")
    expect_equal(actual, expected)
  }
)

context("engroup")

test_that(
  "engroup with capture = TRUE works like capture",
  {
    expected <- capture("foo")
    actual <- engroup("foo", capture = TRUE)
    expect_equal(actual, expected)
  }
)

test_that(
  "engroup with capture = FALSE works like group",
  {
    expected <- group("foo")
    actual <- engroup("foo", capture = FALSE)
    expect_equal(actual, expected)
  }
)

