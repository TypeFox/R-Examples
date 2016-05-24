#options(error = browser)
library(testthat)
library(rpf)

context("sumscore")

test_that("observedSumScore", {
  require(rpf)
  set.seed(1)
  spec <- list()
  spec[1:3] <- rpf.grm(outcomes=3)
  param <- sapply(spec, rpf.rparam, version=1)
  data <- rpf.sample(5, spec, param)
  colnames(param) <- colnames(data)
  grp <- list(spec=spec, param=param, data=data)
  obs <- observedSumScore(grp)
  expect_equal(obs$dist, c(1L, 1L, 0L, 1L, 1L, 0L, 1L))
  
  dperm <- sample.int(3)
  data <- data[,dperm]
  
  mask <- c(TRUE, FALSE, TRUE)
  obs <- observedSumScore(grp, mask=mask)
  expect_equal(obs$dist, rep(1L, 5))
})

test_that("itemOutcomeBySumScore", {
  set.seed(1)
  spec <- list()
  spec[1:3] <- rpf.grm(outcomes=3)
  param <- sapply(spec, rpf.rparam, version=1)
  data <- rpf.sample(5, spec, param)
  colnames(param) <- colnames(data)
  grp <- list(spec=spec, param=param, data=data)
  levels(grp$data[,1]) <- c('a','b','c')
  tbl <- itemOutcomeBySumScore(grp, c(FALSE,TRUE,TRUE), 1L)
  
  want <- structure(c(1L, 0L, 0L, 0L, 0L, 1L, 0L, 0L, 1L, 0L, 0L, 1L, 0L,  0L, 1L),
                    .Dim = c(5L, 3L), .Dimnames = list(c("0", "1", "2",  "3", "4"), c("a","b","c")))
  expect_equal(tbl$table, want)  
})

# Thissen, Pommerich, Billeaud, & Williams (1995)
test_that("tpbw1995-table2", {
  spec <- list()
  spec[1:3] <- rpf.grm(outcomes=4)
  
  param <- matrix(c(1.87, .65, 1.97, 3.14,
                    2.66, .12, 1.57, 2.69,
                    1.24, .08, 2.03, 4.3), nrow=4)
  # fix parameterization
  param <- apply(param, 2, function(p) c(p[1], p[2:4] * -p[1]))
  
  grp <- list(spec=spec, mean=0, cov=matrix(1,1,1), param=param)
  
  got <- sumScoreEAP(grp)
  
  expect_equal(sum(got[,'p']), 1, tolerance=.001)
  
  #cat(deparse(round(got[,2],3)))
  rownames(got) <- NULL
  
  ssP <- c(0.325, 0.241, 0.183, 0.123, 0.069, 0.035, 0.016, 0.006, 0.002,  0)
  expect_equal(got[,'p'], ssP, tolerance=.01)
  ssEAP <- c(-0.885, -0.179, 0.332, 0.744, 1.115, 1.482, 1.843, 2.212, 2.622,  2.999)
  expect_equal(got[,'s1'], ssEAP, tolerance=.01)
  ssVar <- c(0.494, 0.378, 0.329, 0.299, 0.297, 0.296, 0.29, 0.296, 0.313,  0.328)
  expect_equal(got[,'se1'], sqrt(ssVar), tolerance=.01)
  expect_equal(got[,'cov1'], ssVar, tolerance=.01)
})

verifySumP <- function(grp, sseap, N=2000) {  # a good fit is close to 1
  sim <- apply(sapply(rpf.sample(N, grp=grp), unclass), 1, function(r) sum(r-1))
  observed <- tabulate(1+sim, length(sseap[,1]))
  #  print(observed/N)
  ptw2011.gof.test(observed, N*sseap[,1])
}

if (0) {
  fm <- read.flexmirt("~/ifa/ifa-2d-mg/2d-mg-prm.txt")
  
  got <- sumScoreEAP(fm$G1, 5, 21L)  # matches flexmirt exactly
  verifySumP(fm$G1, got)
  
  got <- sumScoreEAP(fm$G2, 5, 21L)  # doesn't match flexmirt
  verifySumP(fm$G2, got, N=5000)  # but looks feasible
  
  got <- sumScoreEAP(fm$G3, 5, 21L)  # doesn't match flexmirt
  verifySumP(fm$G3, got, N=5000)  # but looks feasible
}

test_that("2tier sumScoreEAP", {
  set.seed(1)
  require(rpf)
  numItems <- 6
  spec <- list()
  spec[1:numItems] <- rpf.drm(factors=3)
  param <- sapply(spec, rpf.rparam, version=1)
  gsize <- numItems/3
  for (gx in 0:2) {
    if (gx != 1) {
      param['a2', seq(gx * gsize+1, (gx+1)*gsize)] <- 0
    }
    if (gx != 2) {
      param['a3', seq(gx * gsize+1, (gx+1)*gsize)] <- 0
    }
  }
  grp <- list(spec=spec, param=param, mean=runif(3, -1, 1), cov=diag(runif(3,.5,2)))
  grp$data <- rpf.sample(500, grp=grp)
  colnames(grp$param) <- colnames(grp$data)
  
  got <- sumScoreEAP(grp, qwidth=2, qpoints=5L, .twotier=FALSE)
  tt <- sumScoreEAP(grp, qwidth=2, qpoints=5L, .twotier=TRUE)
  expect_equal(tt, got[,c(1:2,5,8)], .001)
  
  grp2 <- omitItems(grp, c('i5','i6'))
  got <- sumScoreEAP(grp2, qwidth=2, qpoints=5L)
#  cat(deparse(round(log(got[,'p']), 2)))
  expect_equal(log(got[,'p']), c(-3.19, -1.68, -1.07, -1.16, -2.14),
               check.names=FALSE, tolerance=.01)
})
