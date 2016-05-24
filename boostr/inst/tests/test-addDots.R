context("addDots")

test_that("addDots doesn't extend argument when '...' are already there", {
  f <- function(x, y, ...) x + y
  expect_that(addDots(f, .verbose=TRUE), gives_warning())
  expect_that(addDots(f), equals(f))
})

test_that("addDots extends arguments but keeps body intact", {
  f <- function(x, y) x + y
  g <- function(x, y, ...) x + y
  expect_that(addDots(f), equals(g))
  expect_that(body(addDots(f)), equals(body(f)))
})