library(profilr)
context("profile() on vectors")

u <- c(TRUE, TRUE, FALSE, NA)
v <- c("A", "B", "B", NA_character_)
w <- c(0L, 1L, 2L, NA_integer_)
x <- c(5.1, 6.2, 7.3, NA_real_)
y <- seq.Date(from = as.Date("2011-01-01"), to = as.Date("2011-07-01"), length.out = 4)
z <- factor(c("AK", "HI", "CA", "NY"))

test_that("profile() produces data.frames.", {
  expect_is(profile(u), "data.frame")
  expect_is(profile(v), "data.frame")
  expect_is(profile(w), "data.frame")
  expect_is(profile(x), "data.frame")
  expect_is(profile(y), "data.frame")
  expect_is(profile(z), "data.frame")
})
