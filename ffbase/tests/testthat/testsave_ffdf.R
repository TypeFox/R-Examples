library(testthat)

context("save and load")

test_that("Saving and loading works", {
  dir <- tempdir()
  iris_ffdf <- as.ffdf(iris)
  save.ffdf(iris_ffdf, dir=dir)
  rm(iris_ffdf)
  load.ffdf(dir)
})