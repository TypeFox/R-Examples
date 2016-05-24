library(kriens)
context("compose function and forget")

test_that("do(h, g, f) = forget(h %.% g %.% f)", {
  f <- function(x, ret) {
    ret(x+1)
  }
  g <- function(x, ret) {
    ret(x*2)
  }
  h <- function(x, ret){
    ret(x) * ret(x)
  }
  r1 <- forget(h %.% g %.% f)
  r2 <- do(h, g, f)

  for(i in 1:100) {
    expect_equal(r1(i), r2(i))
  }
})
