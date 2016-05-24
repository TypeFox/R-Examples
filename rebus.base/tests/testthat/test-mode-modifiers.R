context("case_insensitive")

test_that(
  "case_insensitive wraps in case insensitivity tokens",
  {
    expected <- as.regex("(?i)foo(?-i)")
    actual <- case_insensitive("foo")
    expect_equal(actual, expected)
  }
)

context("free_spacing")

test_that(
  "free_spacing wraps in free spacing tokens",
  {
    expected <- as.regex("(?x)foo(?-x)")
    actual <- free_spacing("foo")
    expect_equal(actual, expected)
  }
)

context("single_line")

test_that(
  "single_line wraps in single line tokens",
  {
    expected <- as.regex("(?s)foo(?-s)")
    actual <- single_line("foo")
    expect_equal(actual, expected)
  }
)

context("multi_line")

test_that(
  "multi_line wraps in multi-line tokens",
  {
    expected <- as.regex("(?m)foo(?-m)")
    actual <- multi_line("foo")
    expect_equal(actual, expected)
  }
)

context("duplicate_group_names")

test_that(
  "duplicate_group_names wraps in duplicate group name tokens",
  {
    expected <- as.regex("(?J)foo(?-J)")
    actual <- duplicate_group_names("foo")
    expect_equal(actual, expected)
  }
)

context("no_backslash_escaping")

test_that(
  "no_backslash_escaping wraps in no backslash escaping tokens",
  {
    expected <- as.regex("(?X)foo(?-X)")
    actual <- no_backslash_escaping("foo")
    expect_equal(actual, expected)
  }
)
