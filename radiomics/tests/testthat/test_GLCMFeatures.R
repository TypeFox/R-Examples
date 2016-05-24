library(radiomics)
context("GLCM Features")

hbGLCM <- glcm(hallbey, n_grey=4, verbose=FALSE)
test_that("0 degree GLCM Features are properly calculated", {
  #Values taken from http://www.fp.ucalgary.ca/mhallbey/tutorial.htm
  expect_equal(glcm_mean(hbGLCM), 1.291667, tolerance = .02)
  expect_equal(glcm_variance(hbGLCM), 1.039931, tolerance = .02)
  expect_equal(glcm_dissimilarity(hbGLCM), 0.4167, tolerance = .02)
  expect_equal(glcm_contrast(hbGLCM), 0.583, tolerance = .02)
  expect_equal(glcm_homogeneity2(hbGLCM), 0.807, tolerance = 002)
  expect_equal(glcm_energy(hbGLCM), 0.145, tolerance = .02)  
  expect_equal(glcm_entropy(hbGLCM, base=exp(1)), 2.0951, tolerance = .02)
  expect_equal(glcm_correlation(hbGLCM), 0.718, tolerance = .02)  
})

test_that("0 degree GLCM features are properly calculated", {
  expect_equal(calc_features(glcm(hallbey, angle=0, v=F)), read.csv("../csvs/glcm_features/hb0.csv", stringsAsFactors=FALSE))
  expect_equal(calc_features(glcm(tumor, angle=0, v=F)), read.csv("../csvs/glcm_features/tumor0.csv", stringsAsFactors=FALSE))
  expect_equal(calc_features(glcm(noise, angle=0, v=F)), read.csv("../csvs/glcm_features/noise0.csv", stringsAsFactors=FALSE))
  expect_equal(calc_features(glcm(bars, angle=0, v=F)), read.csv("../csvs/glcm_features/bars0.csv", stringsAsFactors=FALSE))
  
})

test_that("45 degree GLCM features are properly calculated", {
  expect_equal(calc_features(glcm(hallbey, angle=45, v=F)), read.csv("../csvs/glcm_features/hb45.csv", stringsAsFactors=FALSE))
  expect_equal(calc_features(glcm(tumor, angle=45, v=F)), read.csv("../csvs/glcm_features/tumor45.csv", stringsAsFactors=FALSE))
  expect_equal(calc_features(glcm(noise, angle=45, v=F)), read.csv("../csvs/glcm_features/noise45.csv", stringsAsFactors=FALSE))
  expect_equal(calc_features(glcm(bars, angle=45, v=F)), read.csv("../csvs/glcm_features/bars45.csv", stringsAsFactors=FALSE))
  
})

test_that("90 degree GLCM features are properly calculated", {
  expect_equal(calc_features(glcm(hallbey, angle=90, v=F)), read.csv("../csvs/glcm_features/hb90.csv", stringsAsFactors=FALSE))
  expect_equal(calc_features(glcm(tumor, angle=90, v=F)), read.csv("../csvs/glcm_features/tumor90.csv", stringsAsFactors=FALSE))
  expect_equal(calc_features(glcm(noise, angle=90, v=F)), read.csv("../csvs/glcm_features/noise90.csv", stringsAsFactors=FALSE))
  expect_equal(calc_features(glcm(bars, angle=90, v=F)), read.csv("../csvs/glcm_features/bars90.csv", stringsAsFactors=FALSE))
  
})

test_that("135 degree GLCM features are properly calculated", {
  expect_equal(calc_features(glcm(hallbey, angle=135, v=F)), read.csv("../csvs/glcm_features/hb135.csv", stringsAsFactors=FALSE))
  expect_equal(calc_features(glcm(tumor, angle=135, v=F)), read.csv("../csvs/glcm_features/tumor135.csv", stringsAsFactors=FALSE))
  expect_equal(calc_features(glcm(noise, angle=135, v=F)), read.csv("../csvs/glcm_features/noise135.csv", stringsAsFactors=FALSE))
  expect_equal(calc_features(glcm(bars, angle=135, v=F)), read.csv("../csvs/glcm_features/bars135.csv", stringsAsFactors=FALSE))
  
})

