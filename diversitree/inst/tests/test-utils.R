source("helper-diversitree.R")

context("Utilities")

test_that("Default argument rewriting works", {
  f <- function(a=1) a
  g <- set.defaults(f, a=2)
  expect_that(formals(f)$a, equals(1))
  expect_that(formals(g)$a, equals(2))
  expect_that(attr(g, "srcref"), is_a("NULL"))
})
