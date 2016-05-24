library(describer)
context("describe() on vectors")

u <- c(TRUE, TRUE, FALSE, NA)
v <- c("A", "B", "B", NA_character_)
w <- c(0L, 1L, 2L, NA_integer_)
x <- c(5.1, 6.2, 7.3, NA_real_)
y <- seq.Date(from = as.Date("2011-01-01"),
              to = as.Date("2011-07-01"),
              length.out = 4)
z <- factor(c("AK", "HI", "CA", "NY"))

test_that("describe() produces data.frames.", {
  expect_is(describe(u), "data.frame")
  expect_is(describe(v), "data.frame")
  expect_is(describe(w), "data.frame")
  expect_is(describe(x), "data.frame")
  expect_is(describe(y), "data.frame")
  expect_is(describe(z), "data.frame")
})
