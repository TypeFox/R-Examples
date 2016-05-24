library(rpf)
library(testthat)

context("null")

test_that("param info", {
  ans1 <-structure(list(type = "intercept", upper = NA_real_, lower = NA_real_),
                   .Names = c("type",  "upper", "lower"))
  expect_identical(rpf.paramInfo(rpf.drm(factors=0)), ans1)
  
  ans2 <- structure(list("intercept", NA_real_, NA_real_, "intercept",      NA_real_, NA_real_),
                    .Dim = c(3L, 2L), .Dimnames = list(c("type",  "upper", "lower"), NULL))
  expect_identical(rpf.paramInfo(rpf.grm(outcomes=3, factors=0)), ans2)
  
  ans3 <- structure(list("intercept", NA_real_, NA_real_,
                         "intercept",      NA_real_, NA_real_,
                         "intercept", NA_real_, NA_real_),
                    .Dim = c(3L,  3L), .Dimnames = list(c("type", "upper", "lower"), NULL))
  expect_identical(rpf.paramInfo(rpf.nrm(outcomes=4, factors=0)), ans3)
})

spec <- list()
param <- list()
spec [[length(spec) +1]] <- rpf.drm(factors = 0)
param[[length(param)+1]] <- c(0)

spec [[length(spec) +1]] <- rpf.grm(outcomes=3, factors=0)
param[[length(param)+1]] <- c(1, -1)

spec [[length(spec) +1]] <- rpf.nrm(outcomes=3, factors=0)
param[[length(param)+1]] <- c(0, -.6)

spec1 <- lapply(spec, rpf.modify, 1)
param1 <- list(c(1,param[[1]],logit(0), logit(1)),
               c(1,param[[2]]),
               c(1,1,0,param[[3]]))

for (ix in 1:length(spec)) {
  test_that(class(spec[[ix]]), {
    expect_equal(rpf.prob(spec[[ix]], param[[ix]], NULL),
                 rpf.prob(spec1[[ix]], param1[[ix]], 0))
    expect_equal(rpf.logprob(spec[[ix]], param[[ix]], NULL),
                 rpf.logprob(spec1[[ix]], param1[[ix]], 0))
    rp <- rpf.rparam(spec[[ix]])
    expect_equal(length(rp), rpf.numParam(spec[[ix]]))
  })
}

test_that("sample null", {
  set.seed(1)
  got <- rpf.sample(3, spec, param)
  got <- sapply(got, unclass)
  colnames(got) <- NULL
  ans <- structure(c(1L, 1L, 1L, 2L, 2L, 2L, 3L, 1L, 3L),
                   .Dim = c(3L,  3L), .Dimnames = list(NULL, NULL))
  expect_equal(got, ans)
})

test_that("null dLL drm", {
  for (rep in 1:3) {
    w <- rchisq(2, 50)
    d0 <- rpf.dLL(spec[[1]], param[[1]], NULL, w)
    d1 <- rpf.dLL(spec1[[1]], param1[[1]], 0, w)
    expect_equal(d0[1], d1[2])
    expect_equal(d0[2], d1[7])
  }
})

test_that("null dLL grm", {
  for (rep in 1:3) {
    w <- rchisq(3, 50)
    item <- 2
    d0 <- rpf.dLL(spec[[item]], param[[item]], NULL, w)
    d1 <- rpf.dLL(spec1[[item]], param1[[item]], 0, w)
    expect_equal(d0[1], d1[2])
    expect_equal(d0[2], d1[3])
    expect_equal(d0[3], d1[6])
    expect_equal(d0[4], d1[8])
    expect_equal(d0[5], d1[9])
  }
})

test_that("null dLL nrm", {
  for (rep in 1:3) {
    w <- rchisq(3, 50)
    item <- 3
    d0 <- rpf.dLL(spec[[item]], param[[item]], NULL, w)
    d1 <- rpf.dLL(spec1[[item]], param1[[item]], 0, w)
    expect_equal(d0[1], d1[4])
    expect_equal(d0[2], d1[5])
    expect_equal(d0[3], d1[15])
    expect_equal(d0[4], 0)
    expect_equal(d0[5], d1[20])
  }
})
