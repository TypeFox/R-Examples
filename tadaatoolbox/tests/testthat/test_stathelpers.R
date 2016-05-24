context("Stat helpers")

df <- data.frame(x = runif(100), y = sample(c("A", "B"), 100, TRUE))

test_that("confint_t returns one numerical value", {
  expect_true(is.numeric(confint_t(df$x)))
  expect_true(length(confint_t(df$x)) == 1)
})

test_that("effect_size returns one numerical value", {
  expect_true(is.numeric(effect_size_t(df, "x", "y")))
  expect_true(length(effect_size_t(df, "x", "y")) == 1)
})
