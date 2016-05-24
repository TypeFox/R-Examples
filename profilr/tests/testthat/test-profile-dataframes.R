library(profilr)
context("profile() on data.frames")

u <- c(TRUE, TRUE, FALSE, NA)
v <- c("A", "B", "B", NA_character_)
w <- c(0L, 1L, 2L, NA_integer_)
x <- c(5.1, 6.2, 7.3, NA_real_)
y <- seq.Date(from = as.Date("2011-01-01"), to = as.Date("2011-07-01"), length.out = 4)
z <- factor(c("AK", "HI", "CA", "NY"))

new_data <- data.frame(u = u, v = v, w = w, x = x, y = y, z = z,
                       stringsAsFactors = FALSE)

test_that("profile() produces data.frames", {
  expect_is(profile(new_data), "data.frame")
})
