source("~/svn/spinyreg/pkg/R/spinyReg.R")

## context("Basic calls to spinyReg function")

## test_that("spinyReg call", {

  ## PROSTATE DATA SET
  load("prostate.rda")
  x <- as.matrix(x)

  out <- spinyreg(x,y,verbose=2)

## })
