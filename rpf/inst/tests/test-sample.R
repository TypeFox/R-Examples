library(testthat)
library(rpf)

context("sample")

unfactor <- function(data) {
  for(s in ls(data)) {
    if (is.factor(data[[s]])) {
      data[[s]] <- unclass(data[[s]])
    }
  }
}

compare.df <- function(df1, df2) {
  df1 <- unfactor(df1)
  df2 <- unfactor(df2)
  expect_identical(dim(df1), dim(df2))
  expect_identical(df1, df2)
}

set.seed(1)

test_that("1 dimensional items", {
  i1 <- rpf.drm()
  i1.p <- rpf.rparam(i1, version=1)
  i2 <- rpf.nrm(outcomes=3)
  i2.p <- rpf.rparam(i2, version=1)
  data <- rpf.sample(3, list(i1,i2), list(i1.p, i2.p))
  compare.df(unclass(data), data.frame(i1=c(2,1,1), i2=c(3,3,3)))

  data <- rpf.sample(runif(3), list(i1,i2), list(i1.p, i2.p))
  compare.df(unclass(data), data.frame(i1=c(2,1,2), i2=c(3,2,3)))
})

test_that("multidimension, no design", {
	set.seed(1)
  numItems <- 3
  items <- vector("list", numItems)
  correct <- vector("list", numItems)

  i1 <- rpf.drm(factors=2)
  for (ix in 1:numItems) {
    items[[ix]] <- i1
    correct[[ix]] <- rpf.rparam(i1, version=1)
  }

  data <- rpf.sample(3, items, correct)
  expect_identical(c(simplify2array(data)),
                   c("2", "1", "2", "2", "1", "1", "1", "2", "2"),
                   info=paste(deparse(simplify2array(data)), collapse="\n"))
})

test_that("1d and 2d", {
  numItems <- 4
  i1 <- rpf.drm()
  i2 <- rpf.drm(factors=2)
  items <- vector("list", numItems)
  for (ix in seq(1,numItems,2)) items[[ix]] <- i1
  for (ix in seq(2,numItems,2)) items[[ix]] <- i2
  correct <- lapply(items, rpf.rparam, version=1)
  data <- rpf.sample(4, items, correct)
  expect_true(all(dim(data) == c(4,4)))
})
