library(kriens)
context("recursively composing functions")

test_that("if an empty list is passed the function errors", {
  expect_error(path())
})

test_that("if NULL is passed the function errors", {
  expect_error(path(NULL))
})

test_that("path(list(h, g, f)) = h %.% g %.% f", {
  f <- function(x, ret) {
    ret(x+1)
  }
  g <- function(x, ret) {
    ret(x*2)
  }
  h <- function(x, ret){
    ret(x) * ret(x)
  }
  r1 <- h %.% g %.% f
  r2 <- path(list(h, g, f))

  for(i in 1:100) {
    expect_equal(r1(i, identity), r2(i, identity))
  }
})
