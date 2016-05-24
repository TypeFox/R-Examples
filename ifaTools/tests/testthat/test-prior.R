library(testthat)
library(ifaTools)

context("prior")

checkAnalyticDeriv <- function(m1) {
  m2 <- mxRun(mxModel(m1, mxComputeSequence(list(
    mxComputeOnce('fitfunction', c('fit','gradient','hessian')),
    mxComputeReportDeriv()))), silent = TRUE)
  m3 <- mxRun(mxModel(m1, mxComputeSequence(list(
    mxComputeNumericDeriv(checkGradient = FALSE),
    mxComputeReportDeriv()))), silent = TRUE)
  
  expect_equal(m2$output$gradient, m3$output$gradient, tolerance=1e-5)
  expect_equal(m2$output$hessian, m3$output$hessian, tolerance=1e-4)
  
  if(0) {
    round(m2$output$hessian,3)
    round(m3$output$hessian,3)
  }
}

test_that("univariatePrior", {
  checkAnalyticDeriv(mxModel(
    univariatePrior("lnorm", "x", 1.1),
    mxMatrix(nrow=1, ncol=1, labels="x", free=TRUE, name="par", values=.3)))

  checkAnalyticDeriv(mxModel(
    univariatePrior("beta", "x", logit(1/4)),
    mxMatrix(nrow=1, ncol=1, labels="x", free=TRUE, name="par", values=logit(.3))))

  checkAnalyticDeriv(mxModel(
    univariatePrior("logit-norm", "x", logit(1/4)),
    mxMatrix(nrow=1, ncol=1, labels="x", free=TRUE, name="par", values=logit(.3))))
})

test_that("uniquenessPrior", {
  imat <- mxMatrix(name='item', nrow=3, ncol=3, values=1:9,
                   free=TRUE, labels=paste("l",1:9,sep=""))
  imat$labels[1,2] <- "l1"
  imat$free[2,2] <- FALSE
  imat$free[3,] <- FALSE
  
  spec <- list()
  spec[1:3] <- rpf.grm(factors=2)
  
  m1 <- mxModel("test", imat,
                mxExpectationBA81(spec))

  m1 <- mxModel(m1, uniquenessPrior(m1, 2),
                mxFitFunctionMultigroup(c("uniquenessPrior")))
  m1$expectation <- NULL

  checkAnalyticDeriv(omxAssignFirstParameters(m1))
})
