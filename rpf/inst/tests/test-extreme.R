library(rpf)
library(testthat)

context("extremes")

test_that("param info", {
  ans1 <- structure(list("slope", NA_real_, 1e-06,
                         "slope", NA_real_, 1e-06,
                         "intercept", NA_real_, NA_real_,
                         "bound", NA_real_, NA_real_,
                         "bound", NA_real_, NA_real_),
                    .Dim = c(3L,  5L), .Dimnames = list(c("type", "upper", "lower"), NULL))
  expect_identical(rpf.paramInfo(rpf.drm(factors=2)), ans1)
  
  ans2 <- structure(list("slope", NA_real_, 1e-06, "slope", NA_real_, 1e-06,
                         "intercept", NA_real_, NA_real_, "intercept", NA_real_, NA_real_),
                    .Dim = 3:4, .Dimnames = list(     c("type", "upper", "lower"), NULL))
  expect_identical(rpf.paramInfo(rpf.grm(outcomes=3, factors=2)), ans2)
  
  ans3 <- structure(list("slope", NA_real_, 1e-06, "slope", NA_real_, 1e-06, "slope",
                         NA_real_, NA_real_, "slope", NA_real_, NA_real_, "slope", NA_real_,
                         NA_real_, "intercept", NA_real_, NA_real_, "intercept", NA_real_, NA_real_,
                         "intercept", NA_real_, NA_real_),
                    .Dim = c(3L, 8L), .Dimnames = list(     c("type", "upper", "lower"), NULL))
  expect_identical(rpf.paramInfo(rpf.nrm(outcomes=4, factors=2)), ans3)
})

spec <- list()
param <- list()
# repair the poor version of drm TODO
#spec [[length(spec) +1]] <- rpf.drm(poor=TRUE)
#param[[length(param)+1]] <- c(1, 0, 0)
spec [[length(spec) +1]] <- rpf.drm()
param[[length(param)+1]] <- c(1, 0, logit(.05), logit(.95))

spec [[length(spec) +1]] <- rpf.grm(3)
param[[length(param)+1]] <- c(1, 1, -1)

spec [[length(spec) +1]] <- rpf.nrm(3)
param[[length(param)+1]] <- c(1,  .5, .6, 0, -.6)

# To debug, set breakpoint on Rf_error 
where <- seq(-1000, 1000, 100)
for (ix in 1:length(spec)) {
  ispec <- spec[[ix]]
  iparam <- param[[ix]]
  test_that(paste("extreme values in", class(ispec)), {
    for (wh in where) {
      v <- rpf.prob(ispec, iparam, wh)
      expect_equal(sum(v), 1)
      if (wh != 0) {
        sep <- sort(-v)[1:2]
        expect_true(abs(sep[1] - sep[2]) > .89,
                    info=paste(c(wh, sep), collapse=" "))
      }
    }
    for (wh in where) {
      v <- rpf.logprob(ispec, iparam, wh)
      if (all(v > -35 & v < 35)) {
        expect_equal(sum(exp(v)), 1)
      }
    }
    w <- rchisq(ispec$outcomes, df=6)
    for (wh in where) rpf.dLL(ispec, iparam, wh, w)
  })
}

for (ix in 1:length(spec)) {
  ispec <- spec[[ix]]
  iparam <- param[[ix]]
  test_that(paste("score=NA", class(ispec)), {
    v <- rpf.prob(ispec, iparam, as.numeric(c(NA,NA)))
    expect_true(all(is.na(v)))
    v <- rpf.logprob(ispec, iparam, as.numeric(c(NA,NA)))
    expect_true(all(is.na(v)))
  })
}
