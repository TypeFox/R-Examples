context("Rate curves")

test_that("One rate only", {  
  r <- 0.1
  d <- zero_eff_to_disc(r)
  expect_that(d, equals(1 / (1 + r)))
  expect_that(disc_to_zero_eff(d), equals(r))
  expect_that(disc_to_french(d), equals(r))
  expect_that(disc_to_german(d), equals(r))
  expect_that(disc_to_fut(d), equals(r))
  expect_that(disc_to_swap(d), equals(r))
})

test_that("Flat curve", {
  f <- rep.int(x = 0.1, times = 12)
  d <- fut_to_disc(f)
  expect_that(d, equals(1 / cumprod(1 + f)))
  expect_that(disc_to_zero_eff(d), equals(f))
  expect_that(disc_to_french(d), equals(f))
  expect_that(disc_to_german(d), equals(f))
  expect_that(disc_to_fut(d), equals(f))
  expect_that(disc_to_swap(d), equals(f))
})

test_that("Upward curve", {
  f <- c(0.1, 0.2, 0.3)
  d <- fut_to_disc(f)
  expect_that(d, equals(1 / cumprod(1 + f)))
  expect_that(disc_to_fut(d), equals(f))
  expect_that(disc_to_zero_eff(d), equals(cumprod(1 + f) ^ (1 / seq_along(f)) - 1))
  rates <- vapply(
    X = seq_along(f), 
    FUN = function(m) find_rate(m = m, d = d, loan_type = "french"), 
    FUN.VALUE = 1)
  expect_that(disc_to_french(d), equals(rates))
  rates <- vapply(
    X = seq_along(f), 
    FUN = function(m) find_rate(m = m, d = d, loan_type = "german"), 
    FUN.VALUE = 1)
  expect_that(disc_to_german(d), equals(rates))
  rates <- vapply(
    X = seq_along(f), 
    FUN = function(m) find_rate(m = m, d = d, loan_type = "bullet"), 
    FUN.VALUE = 1)
  expect_that(disc_to_swap(d), equals(rates))  
})