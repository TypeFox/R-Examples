library(kriens)
context("function composition")

test_that("composition rule: f %.% identity2 = f", {
  f <- function(x, ret) {
    ret(x+1)*2
  }
  g <- f %.% identity2

  for(i in 1:100) {
    expect_equal(g(i,identity), f(i,identity))
  }
})

test_that("composition rule: h %.% (g %.% f) = (h %.% g) %.% f", {
  f <- function(x, ret) {
    ret(x+1)
  }
  g <- function(x, ret) {
    ret(x*2)
  }
  h <- function(x, ret){
    ret(x) * ret(x)
  }
  r1 <- h %.% (g %.% f)
  r2 <- (h %.% g) %.% f

  for(i in 1:100) {
    expect_equal(r1(i, identity), r2(i, identity))
  }
})
