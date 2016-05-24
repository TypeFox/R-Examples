################################################################################
# Internal functions
################################################################################

context("internal functions")

test_that("trim_ws works properly", {
  x <- c("abc", " abc", "abc ", " abc ")
  y <- c("123", " 123", "123 ", " 123 ")
  d <- data.frame(x, y)
  trim_both <- trim_ws(x, which = "both")
  trim_left <- trim_ws(x, which = "left")
  trim_right <- trim_ws(x, which = "right")
  expect_identical(trim_both, c("abc", "abc", "abc", "abc"))
  expect_identical(trim_left, c("abc", "abc", "abc ", "abc "))
  expect_identical(trim_right, c("abc", " abc", "abc", " abc"))
  expect_equal(trim_ws(d), data.frame(x = rep("abc", 4), y = rep("123", 4),
                                      stringsAsFactors = FALSE))
})

