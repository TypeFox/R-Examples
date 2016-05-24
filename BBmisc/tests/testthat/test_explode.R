context("explode")

test_that("explode", {
  x = "R is a nice programming language"
  substrings = c("R", "is", "a", "nice", "programming", "language")
  sep = " "

  # split string
  exploded = explode(x, sep = sep)
  expect_equal(length(exploded), 6)
  for (i in 1:length(substrings)) {
    expect_equal(substrings[i], exploded[[i]])
  }

  # now glue the substrings together
  collapsed = collapse(exploded, sep = sep)
  expect_equal(collapsed, x)
})
