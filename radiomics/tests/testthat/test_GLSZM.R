library(radiomics)
context("GLSZM")

test_that("GLSZM is properly generated", {
  expect_error(glszm(data.frame()))
  expect_error(glszm())
  expect_error(glszm(c()))
  expect_error(glszm(list()))
  expect_warning(glszm(matrix()))
  expect_is(glszm(hallbey, n_grey=4), "glszm")
  expect_is(glszm(hallbey, n_grey=4), "matrix")
})

test_that("GLSZM properly calculates", {
  expect_equal(glszm(hallbey, v=F)@.Data, as.matrix(read.table("../csvs/glszm/hb.csv", header=TRUE, sep=",", check.names=FALSE)))
  expect_equal(glszm(tumor, v=F)@.Data, as.matrix(read.table("../csvs/glszm/tumor.csv", header=TRUE, sep=",", check.names=FALSE)))
  expect_equal(glszm(noise, v=F)@.Data, as.matrix(read.table("../csvs/glszm/noise.csv", header=TRUE, sep=",", check.names=FALSE)))
  expect_equal(glszm(bars, v=F)@.Data, as.matrix(read.table("../csvs/glszm/bars.csv", header=TRUE, sep=",", check.names=FALSE)))
  
})