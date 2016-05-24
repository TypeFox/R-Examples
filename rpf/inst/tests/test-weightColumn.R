#options(error = browser)
library(testthat)
library(rpf)

context("weightColumn")

data(LSAT6)
LSAT6[1:5] <- as.data.frame(lapply(LSAT6[1:5], ordered, 0:1))
LSAT6$Freq <- as.numeric(LSAT6$Freq)
spec <- list()
spec[1:5] <- rpf.grm()
names(spec) <- colnames(LSAT6)[1:5]
param <- matrix(c(1.0,0.0), nrow=2, ncol=5, dimnames=list(c('a','b'), names(spec)))

grp1 <- list(spec=spec,
             param=param,
             data=LSAT6,
             weightColumn='Freq',
             observedStats=nrow(LSAT6)-1L,
             minItemsPerScore=5L)
grp1$score <- EAPscores(grp1)

grp2 <- grp1
grp2$data <- expandDataFrame(LSAT6, freqName = "Freq")
grp2$weightColumn <- NULL
grp2$score <- EAPscores(grp2)

expect_equal(nrow(EAPscores(grp1, compressed=TRUE)), nrow(LSAT6))
expect_equal(nrow(grp1$score), nrow(grp2$score))
expect_equal(grp1$score[,1], grp2$score[,1])

test_that("observedSumScore", {
  t1 <- observedSumScore(grp1)
  t2 <- observedSumScore(grp2)
  expect_equal(t1$n, t2$n)
  expect_equal(t1$dist, t2$dist)
})

test_that("sumScoreEAPTest", {
  t1 <- sumScoreEAPTest(grp1)
  t2 <- sumScoreEAPTest(grp2)
  expect_equal(t1$n, t2$n)
  expect_equal(t1$observed, t2$observed)
  expect_equal(t1$pearson.p, t2$pearson.p)
})

test_that("itemOutcomeBySumScore", {
  t1 <- itemOutcomeBySumScore(grp1, c(FALSE, rep(TRUE,4)), 1)
  t2 <- itemOutcomeBySumScore(grp2, c(FALSE, rep(TRUE,4)), 1)
  expect_equal(t1$n, t2$n)
  expect_equal(t1$table, t2$table)
})

test_that("ChenThissen1997", {
  t1 <- ChenThissen1997(grp1)
  t2 <- ChenThissen1997(grp2)
  expect_equal(t1$pval, t2$pval)
})

test_that("SitemFit", {
  t1 <- SitemFit(grp1)
  tbl1 <- t(sapply(t1, function(r) c(n=r$n, df=r$df, stat=r$statistic, pval=r$pval)))
  t2 <- SitemFit(grp2)
  tbl2 <- t(sapply(t2, function(r) c(n=r$n, df=r$df, stat=r$statistic, pval=r$pval)))
  expect_equal(tbl1, tbl2)
})
