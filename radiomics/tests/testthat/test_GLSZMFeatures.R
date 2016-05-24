library(radiomics)
context("GLSZM Features")

hbGLSZM <- glszm(hallbey, n_grey=4)
test_that("GLSZM Features are properly calculated", {
  #Values taken from http://www.fp.ucalgary.ca/mhallbey/tutorial.htm
  expect_equal(glszm_SAE(hbGLSZM), 0.098125, tolerance = .02)
  expect_equal(glszm_LAE(hbGLSZM), 17.5, tolerance = .02)
  expect_equal(glszm_IV(hbGLSZM), 1, tolerance = .02)
  expect_equal(glszm_SZV(hbGLSZM), 1.5, tolerance = .02)
  expect_equal(glszm_ZP(hbGLSZM), 0.25, tolerance = 002)
  expect_equal(glszm_LIE(hbGLSZM), 0.34, tolerance = .02)  
  expect_equal(glszm_HIE(hbGLSZM), 3.5, tolerance = .02)
  expect_equal(glszm_LISAE(hbGLSZM), 0.025, tolerance = .02)  
  expect_equal(glszm_HISAE(hbGLSZM), 0.618, tolerance = .02)
  expect_equal(glszm_LILAE(hbGLSZM), 5.67, tolerance = .02) 
  expect_equal(glszm_HILAE(hbGLSZM), 38, tolerance = .02)
})

test_that("GLSZM features are properly calculated", {
  expect_equal(calc_features(glszm(hallbey, v=F)), read.csv("../csvs/glszm_features/hb.csv", stringsAsFactors=FALSE))
  expect_equal(calc_features(glszm(tumor, v=F)), read.csv("../csvs/glszm_features/tumor.csv", stringsAsFactors=FALSE))
  expect_equal(calc_features(glszm(noise, v=F)), read.csv("../csvs/glszm_features/noise.csv", stringsAsFactors=FALSE))
  expect_equal(calc_features(glszm(bars, v=F)), read.csv("../csvs/glszm_features/bars.csv", stringsAsFactors=FALSE))
  
})