context("or, not enough args")

test_that(
  "or with zero args in ... gives a warning",
  {
    expected <- as.regex("(?:)")
    expect_warning(
      actual <- or(),
      "'or' is intended to be called with at least 2 arguments in '...'. 0 were passed. Maybe you wanted 'or1' instead?"
    )
    expect_equal(actual, expected)
  }
)

test_that(
  "or with one arg in ... gives a warning",
  {
    expected <- as.regex("(?:foo)")
    expect_warning(
      actual <- or("foo"),
      "'or' is intended to be called with at least 2 arguments in '...'. 1 was passed. Maybe you wanted 'or1' instead?"
    )
    expect_equal(actual, expected)
  }
)

test_that(
  "or with two args including capture gives a warning",
  {
    expected <- as.regex("(foo)")
    expect_warning(
      actual <- or("foo", capture = TRUE),
      "'or' is intended to be called with at least 2 arguments in '...'. 1 was passed. Maybe you wanted 'or1' instead?"
    )
    expect_equal(actual, expected)
  }
)

context("or, capture arg")

test_that(
  "or with >= two args separates with |, groups",
  {
    expected <- as.regex("(?:foo|bar|baz)")
    actual <- or("foo", "bar", "baz")
    expect_equal(actual, expected)
  }
)

test_that(
  "or with capture = TRUE captures",
  {
    expected <- as.regex("(foo|bar|baz)")
    actual <- or("foo", "bar", "baz", capture = TRUE)
    expect_equal(actual, expected)
  }
)

test_that(
  "or with capture = NA doesn't group",
  {
    expected <- as.regex("foo|bar|baz")
    actual <- or("foo", "bar", "baz", capture = NA)
    expect_equal(actual, expected)
  }
)

context("or1")

test_that(
  "or1 with capture = TRUE captures",
  {
    expected <- as.regex("(foo|bar|baz)")
    actual <- or1(c("foo", "bar", "baz"), capture = TRUE)
    expect_equal(actual, expected)
  }
)

test_that(
  "or1 with capture = NA doesn't group",
  {
    expected <- as.regex("foo|bar|baz")
    actual <- or1(c("foo", "bar", "baz"), capture = NA)
    expect_equal(actual, expected)
  }
)

context("%|%")

test_that(
  "%|% separates with |, doesn't group",
  {
    expected <- as.regex("foo|bar")
    actual <- "foo" %|% "bar"
    expect_equal(actual, expected)
  }
)
