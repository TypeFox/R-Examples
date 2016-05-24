context("%R%")

test_that(
  "%R% concatenates",
  {
    expected <- as.regex("foobar")
    actual <- "foo" %R% "bar"
    expect_equal(actual, expected)
  }
)
