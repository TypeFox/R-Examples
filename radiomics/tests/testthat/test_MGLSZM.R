library(radiomics)
context("MGLSZM")

test_that("MGLSZM is properly generated", {
  expect_error(mglszm(data.frame()))
  expect_error(mglszm())
  expect_error(mglszm(c()))
  expect_error(mglszm(list()))
  expect_warning(mglszm(matrix()))
  expect_is(mglszm(hallbey), "mglszm")
  expect_is(mglszm(hallbey), "matrix")
})

test_that("MGLSZM properly calculates", {
  expect_equal(mglszm(hallbey)@.Data, as.matrix(read.table("../csvs/mglszm/hb.csv", header=TRUE, sep=",", check.names=FALSE)))
  expect_equal(mglszm(tumor)@.Data, as.matrix(read.table("../csvs/mglszm/tumor.csv", header=TRUE, sep=",", check.names=FALSE)))
  expect_equal(mglszm(noise)@.Data, as.matrix(read.table("../csvs/mglszm/noise.csv", header=TRUE, sep=",", check.names=FALSE)))
  expect_equal(mglszm(bars)@.Data, as.matrix(read.table("../csvs/mglszm/bars.csv", header=TRUE, sep=",", check.names=FALSE)))
  
})