library(testthat)
library(Causata)

equals <- testthat::equals

context("GetStratifiedSample")

test_that("basic sampling probabilities are correct with more negative", {
    df <- data.frame(iv=rnorm(10), dv=c(0, 1, 0, 1, 1, 1, 0, 0, 0, 0))
    sample.probabilities <- Causata:::SamplingProbabilities(df$dv, 0)
    
    expect_that(length(sample.probabilities), equals(2))
    expect_that(sample.probabilities$A, equals(0.8164966, tolerance=0.000001))
    expect_that(sample.probabilities$B, equals(1))
  } 
)

test_that("basic sampling probabilities are correct with more positive", {
    df <- data.frame(iv=rnorm(10), dv=c(1, 1, 0, 1, 1, 0, 1, 0,1, 1))
    sample.probabilities <- Causata:::SamplingProbabilities(df$dv, 0)
    
    expect_that(length(sample.probabilities), equals(2))
    expect_that(sample.probabilities$A, equals(1))
    expect_that(sample.probabilities$B, equals(0.6546537, tolerance=0.000001))
  } 
)
          
test_that("sampling where stratification.variable is single-valued and negative", {
    df <- data.frame(iv=rnorm(10), dv=c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    sample.probabilities <- Causata:::SamplingProbabilities(df$dv, 0)
    
   expect_that(length(sample.probabilities), equals(2))
   expect_that(sample.probabilities$A, equals(0.3162278, tolerance=0.000001))
   expect_that(sample.probabilities$B, equals(1))
  }
)

test_that("sampling where stratification.variable is single-valued and positive", {
    df <- data.frame(iv=rnorm(10), dv=c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1))
    sample.probabilities <- Causata:::SamplingProbabilities(df$dv, 0)
    
    expect_that(length(sample.probabilities), equals(2))
    expect_that(sample.probabilities$A, equals(1))
    expect_that(sample.probabilities$B, equals(0.3162278, tolerance=0.000001))
  }
)

test_that("stratification value", {
    df <- data.frame(iv=rnorm(10), dv=c(0, 1, 0, 3, 1, 1, 0, 0, 0, 3))
    sample.probabilities <- Causata:::SamplingProbabilities(df$dv, stratification.value=3)
    
    expect_that(length(sample.probabilities), equals(2))
    expect_that(sample.probabilities$A, equals(1))
    expect_that(sample.probabilities$B, equals(0.5))
  }
)

test_that("stratification value is floating point", {
    df <- data.frame(iv=rnorm(10), dv=c(0, 0.012, 0, 0.234324, 134.34, 13.4, 0, 0, 0, 0))
    sample.probabilities <- Causata:::SamplingProbabilities(df$dv, 0)
    
    expect_that(length(sample.probabilities), equals(2))
    expect_that(sample.probabilities$A, equals(0.8164966, tolerance=0.000001))
    expect_that(sample.probabilities$B, equals(1))
  }
)