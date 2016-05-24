library(radiomics)
context("GLCM")

test_that("GLCM is handles different inputs properly", {
  expect_error(glcm(data.frame()))
  expect_error(glcm())
  expect_error(glcm(c()))
  expect_error(glcm(list()))
  expect_warning(glcm(matrix()))
  expect_is(glcm(hallbey, n_grey=4), "glcm")
  expect_is(glcm(hallbey, n_grey=4), "matrix")
})

test_that("0 degree GLCM properly calculates", {
  expect_equal(glcm(hallbey, angle=0, v=F)@.Data, as.matrix(read.table("../csvs/glcm/hb0.csv", header=TRUE, sep=",", check.names=FALSE)))
  expect_equal(glcm(tumor, angle=0, v=F)@.Data, as.matrix(read.table("../csvs/glcm/tumor0.csv", header=TRUE, sep=",", check.names=FALSE)))
  expect_equal(glcm(noise, angle=0, v=F)@.Data, as.matrix(read.table("../csvs/glcm/noise0.csv", header=TRUE, sep=",", check.names=FALSE)))
  expect_equal(glcm(bars, angle=0, v=F)@.Data, as.matrix(read.table("../csvs/glcm/bars0.csv", header=TRUE, sep=",", check.names=FALSE)))
  
})

test_that("45 degree GLCM properly calculates", {
  expect_equal(glcm(hallbey, angle=45, v=F)@.Data, as.matrix(read.table("../csvs/glcm/hb45.csv", header=TRUE, sep=",", check.names=FALSE)))
  expect_equal(glcm(tumor, angle=45, v=F)@.Data, as.matrix(read.table("../csvs/glcm/tumor45.csv", header=TRUE, sep=",", check.names=FALSE)))
  expect_equal(glcm(noise, angle=45, v=F)@.Data, as.matrix(read.table("../csvs/glcm/noise45.csv", header=TRUE, sep=",", check.names=FALSE)))
  expect_equal(glcm(bars, angle=45, v=F)@.Data, as.matrix(read.table("../csvs/glcm/bars45.csv", header=TRUE, sep=",", check.names=FALSE)))
  
})

test_that("90 degree GLCM properly calculates", {
  expect_equal(glcm(hallbey, angle=90, v=F)@.Data, as.matrix(read.table("../csvs/glcm/hb90.csv", header=TRUE, sep=",", check.names=FALSE)))
  expect_equal(glcm(tumor, angle=90, v=F)@.Data, as.matrix(read.table("../csvs/glcm/tumor90.csv", header=TRUE, sep=",", check.names=FALSE)))
  expect_equal(glcm(noise, angle=90, v=F)@.Data, as.matrix(read.table("../csvs/glcm/noise90.csv", header=TRUE, sep=",", check.names=FALSE)))
  expect_equal(glcm(bars, angle=90, v=F)@.Data, as.matrix(read.table("../csvs/glcm/bars90.csv", header=TRUE, sep=",", check.names=FALSE)))
  
})

test_that("135 degree GLCM properly calculates", {
  expect_equal(glcm(hallbey, angle=135, v=F)@.Data, as.matrix(read.table("../csvs/glcm/hb135.csv", header=TRUE, sep=",", check.names=FALSE)))
  expect_equal(glcm(tumor, angle=135, v=F)@.Data, as.matrix(read.table("../csvs/glcm/tumor135.csv", header=TRUE, sep=",", check.names=FALSE)))
  expect_equal(glcm(noise, angle=135, v=F)@.Data, as.matrix(read.table("../csvs/glcm/noise135.csv", header=TRUE, sep=",", check.names=FALSE)))
  expect_equal(glcm(bars, angle=135, v=F)@.Data, as.matrix(read.table("../csvs/glcm/bars135.csv", header=TRUE, sep=",", check.names=FALSE)))
  
})
