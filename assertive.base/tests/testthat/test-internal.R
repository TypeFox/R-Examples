test_that(
  "print_and_capture returns a string of the captured output.",
  {
    expected <- "  x\n1 a\n2 b\n3 c"
    actual <- print_and_capture(data.frame(x = letters[1:3]))
    expect_equal(actual, expected)
  }
)
