context("Check inv.logit and logit")

p <- runif(100)
a <- GMCM:::logit(p)
pp <- GMCM:::inv.logit(a)

test_that("inv.logit and logit returns proper formatted output", {
  expect_that(is.numeric(a),  is_true())
  expect_that(length(a),  equals(length(p)))
  expect_that(is.numeric(pp),  is_true())
  expect_that(length(pp),  equals(length(p)))
  expect_that(p - pp,  equals(rep(0, length(p))))
})

# Test degenerate input
