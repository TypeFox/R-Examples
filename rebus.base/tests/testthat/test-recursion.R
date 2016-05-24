context("recursive")

test_that(
  "recursive append recursion token",
  {
    expected <- as.regex("foo(?R)")
    actual <- recursive("foo")
    expect_equal(actual, expected)
  }
)

