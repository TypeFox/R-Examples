#options(error = browser)
library(testthat)
library(rpf)

context("score")

test_that("tpbw1995-table2", {
  set.seed(1)
  spec <- list()
  spec[1:3] <- rpf.grm(outcomes=4)
  
  param <- matrix(c(1.87, .65, 1.97, 3.14,
                    2.66, .12, 1.57, 2.69,
                    1.24, .08, 2.03, 4.3), nrow=4)
  # fix parameterization
  param <- apply(param, 2, function(p) c(p[1], p[2:4] * -p[1]))
  colnames(param) <- paste('i', 1:3, sep="")
  
  grp <- list(spec=spec, mean=0, cov=matrix(1,1,1), param=param)
  
  expect_error(EAPscores(grp), "EAP requested but there are no data rows")

  grp$data <- rpf.sample(2, spec, param)
  scores <- EAPscores(grp)
  
  expect_equal(scores[,1], c(0.084, -0.154), tolerance=.01)
  expect_equal(scores[,2], c(0.539, 0.572), tolerance=.01)
})
