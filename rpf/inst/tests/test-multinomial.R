#options(error = browser)
library(testthat)
library(rpf)

context("multinomialFit")

test_that("multinomialFit full info, simple", {
  require(rpf)
  set.seed(7)
  grp <- list(spec=list())
  grp$spec[1:10] <- rpf.grm()
  grp$param <- sapply(grp$spec, rpf.rparam)
  colnames(grp$param) <- paste("i", 1:10, sep="")
  grp$mean <- 0
  grp$cov <- diag(1)
  grp$uniqueFree <- sum(grp$param != 0)
  grp$data <- compressDataFrame(rpf.sample(1500, grp=grp))
  grp$weightColumn <- 'freq'
  grp$observedStats <- nrow(grp$data) - 1
  
  # screw up fit
  grp$param[2,] <- grp$param[2,] + runif(10, -.75, .75)

  got <- multinomialFit(grp, method="pearson")
  expect_equal(got$statistic, 1287.46, tolerance=.01)
  expect_equal(got$df, 209)
  expect_equal(got$pval, -354, tolerance=.1)

  got <- multinomialFit(grp, method="lr")
  expect_equal(got$statistic, 914.35, tolerance=.01)
  expect_equal(got$df, 209)
  expect_equal(got$pval, -202, tolerance=.1)
})

test_that("multinomialFit full info, simple w/ missingness", {
  require(rpf)
  require(testthat)
  
  set.seed(7)
  grp <- list(spec=list())
  grp$spec[1:10] <- rpf.grm()
  grp$param <- sapply(grp$spec, rpf.rparam)
  colnames(grp$param) <- paste("i", 1:10, sep="")
  grp$mean <- 0
  grp$cov <- diag(1)
  grp$free <- grp$param != 0
  grp$labels <- matrix(NA, nrow(grp$param), ncol(grp$param))
  grp$uniqueFree <- sum(grp$param != 0)
  grp$data <- compressDataFrame(rpf.sample(1500, grp=grp, mcar=.1))
  grp$weightColumn <- 'freq'
  grp$observedStats <- nrow(grp$data) - 1
  
  got <- multinomialFit(grp, method="pearson")
  expect_equal(got$statistic, 1115.63, tolerance=.01)
  expect_equal(got$df, 710)
  expect_equal(got$pval, -45.7, tolerance=.1)

  got <- multinomialFit(grp, method="lr")
  expect_equal(got$statistic, 414.691, tolerance=.01)
  expect_equal(got$df, 710)
  expect_equal(got$pval, 0, tolerance=.01)
  expect_equal(got$n, 521)

  got <- multinomialFit(omitMostMissing(grp, 1), method="lr")
  expect_equal(got$statistic, 253.39, tolerance=.01)
  expect_equal(got$df, 585)
  expect_equal(got$pval, 0, tolerance=.01)
  expect_equal(got$n, 598)
})

test_that("multinomialFit full info, two-tier", {
  require(rpf)
  spec <- list()
  spec[1:5] <- rpf.drm(factors=3)
  gen.param <- sapply(spec, rpf.rparam)
  gen.param['a2', 1:2] <- 0
  gen.param['a3', 3] <- 0
  gen.param[c('a2','a3'), 4:5] <- 0
  colnames(gen.param) <- paste("i", 1:ncol(gen.param), sep="")
  resp <- rpf.sample(1500, spec, gen.param)
  grp <- list(spec=spec, param=gen.param, mean=runif(3, 0, 1), cov=diag(runif(3,1,2)),
              data=resp, uniqueFree=sum(gen.param!=0), qwidth=5, qpoints=21)
  
  got1 <- multinomialFit(grp, .twotier = FALSE)
  got2 <- multinomialFit(grp, .twotier = TRUE)
  expect_equal(got1$statistic, got2$statistic, tolerance=.001)
  expect_equal(got1$df, got2$df)
  expect_equal(got1$pval, got2$pval, tolerance=.001)
})

if (0) {
  # matches flexmirt exactly
  library(mirt)
  dat <- expand.table(LSAT7)
  (mod1 <- mirt(dat, 1, calcNull=TRUE))
  i1 <- extract.item(mod1, 1)
  expected.item(i1, seq(-4,4,.5), min=0L)
  M2(mod1)
  M2(mod1@null.mod)
  write.table(dat, file = "LSAT7.csv", row.names = FALSE, col.names=FALSE)
}
