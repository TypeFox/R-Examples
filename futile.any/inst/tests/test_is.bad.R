context("is.bad")

test_that("list with bad values", {
  a <- list(a=1:3, b=NULL, c=NA, d='foo')
  e <- list(a=rep(FALSE,3), b=TRUE, c=TRUE, d=FALSE)
  expect_that(is.bad(a), equals(e) )
})

test_that("vector with NAs", {
  a <- c(1,NA,3)
  expect_that(is.bad(a), equals(c(FALSE,TRUE,FALSE)))
})

test_that("data.frame with NAs", {
  a <- data.frame(a=1:3, b=NA)
  e <- matrix(c(rep(FALSE,3), rep(TRUE,3)), ncol=2)
  colnames(e) <- c('a','b')
  expect_that(is.bad(a), equals(e))
})

test_that("data.frame that is empty", {
  a <- data.frame(a=NULL, b=NULL)
  expect_that(is.bad(a), equals(TRUE))
})

test_that("matrix with NAs", {
  a <- matrix(c(1:3, NA), ncol=2)
  e <- matrix(c(rep(FALSE,3), TRUE), ncol=2)
  expect_that(is.bad(a), equals(e))
})

