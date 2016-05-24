context("escape_special")

test_that(
  "escape_special escapes all special values",
  {
    expected <- as.regex("\\\\ \\^ \\$ \\. \\| \\? \\* \\+ \\( \\) \\{ } \\[ ]")
    actual <- escape_special("\\ ^ $ . | ? * + ( ) { } [ ]")
    expect_equal(actual, expected)
  }
)
