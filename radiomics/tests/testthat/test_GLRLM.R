library(radiomics)
context("GLRLM")

test_that("GLRLM is properly generated", {
  expect_error(glrlm(data.frame()))
  expect_error(glrlm())
  expect_error(glrlm(c()))
  expect_error(glrlm(list()))
  expect_warning(glrlm(matrix()))
  expect_is(glrlm(hallbey, n_grey=4), "glrlm")
  expect_is(glrlm(hallbey, n_grey=4), "matrix")
})

test_that("0 degree GLRLM properly calculates", {
  expect_equal(glrlm(hallbey, angle=0, v=F)@.Data, as.matrix(read.table("../csvs/glrlm/hb0.csv", header=TRUE, sep=",", check.names=FALSE)))
  expect_equal(glrlm(tumor, angle=0, v=F)@.Data, as.matrix(read.table("../csvs/glrlm/tumor0.csv", header=TRUE, sep=",", check.names=FALSE)))
  expect_equal(glrlm(noise, angle=0, v=F)@.Data, as.matrix(read.table("../csvs/glrlm/noise0.csv", header=TRUE, sep=",", check.names=FALSE)))
  expect_equal(glrlm(bars, angle=0, v=F)@.Data, as.matrix(read.table("../csvs/glrlm/bars0.csv", header=TRUE, sep=",", check.names=FALSE)))
  
})

test_that("45 degree GLRLM properly calculates", {
  expect_equal(glrlm(hallbey, angle=45, v=F)@.Data, as.matrix(read.table("../csvs/glrlm/hb45.csv", header=TRUE, sep=",", check.names=FALSE)))
  expect_equal(glrlm(tumor, angle=45, v=F)@.Data, as.matrix(read.table("../csvs/glrlm/tumor45.csv", header=TRUE, sep=",", check.names=FALSE)))
  expect_equal(glrlm(noise, angle=45, v=F)@.Data, as.matrix(read.table("../csvs/glrlm/noise45.csv", header=TRUE, sep=",", check.names=FALSE)))
  expect_equal(glrlm(bars, angle=45, v=F)@.Data, as.matrix(read.table("../csvs/glrlm/bars45.csv", header=TRUE, sep=",", check.names=FALSE)))
  
})

test_that("0 degree GLRLM properly calculates", {
  expect_equal(glrlm(hallbey, angle=90, v=F)@.Data, as.matrix(read.table("../csvs/glrlm/hb90.csv", header=TRUE, sep=",", check.names=FALSE)))
  expect_equal(glrlm(tumor, angle=90, v=F)@.Data, as.matrix(read.table("../csvs/glrlm/tumor90.csv", header=TRUE, sep=",", check.names=FALSE)))
  expect_equal(glrlm(noise, angle=90, v=F)@.Data, as.matrix(read.table("../csvs/glrlm/noise90.csv", header=TRUE, sep=",", check.names=FALSE)))
  expect_equal(glrlm(bars, angle=90, v=F)@.Data, as.matrix(read.table("../csvs/glrlm/bars90.csv", header=TRUE, sep=",", check.names=FALSE)))
  
})

test_that("45 degree GLRLM properly calculates", {
  expect_equal(glrlm(hallbey, angle=135, v=F)@.Data, as.matrix(read.table("../csvs/glrlm/hb135.csv", header=TRUE, sep=",", check.names=FALSE)))
  expect_equal(glrlm(tumor, angle=135, v=F)@.Data, as.matrix(read.table("../csvs/glrlm/tumor135.csv", header=TRUE, sep=",", check.names=FALSE)))
  expect_equal(glrlm(noise, angle=135, v=F)@.Data, as.matrix(read.table("../csvs/glrlm/noise135.csv", header=TRUE, sep=",", check.names=FALSE)))
  expect_equal(glrlm(bars, angle=135, v=F)@.Data, as.matrix(read.table("../csvs/glrlm/bars135.csv", header=TRUE, sep=",", check.names=FALSE)))
  
})