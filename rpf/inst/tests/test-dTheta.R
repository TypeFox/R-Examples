library(testthat)
library(rpf)
library(numDeriv)
#options(error = utils::recover)

context("dTheta")

items <- list(rpf.drm(),
              rpf.drm(factors=2),
              rpf.grm(outcomes=3, factors=2),
              rpf.nrm(outcomes=4, factors=2, T.a="random", T.c="random"),
              rpf.nrm(outcomes=4, factors=3, T.a="random", T.c="random"),
              rpf.lmp(k=0),
              rpf.lmp(k=1))

triSize <- function(sz) sz*(sz + 1) / 2

for (ii in items) {
  test_that(class(ii), {
    dir <- runif(ii$factors)
    dir <- dir / sqrt(sum(dir^2))

    ii.p <- rpf.rparam(ii)
    at <- rnorm(ii$factors)
    analytic <- rpf.dTheta(ii, ii.p, at, dir)
    numDeriv <- ii$factors + triSize(ii$factors)
    deriv <- matrix(NA, numDeriv, ii$outcomes)
    for (ox in 1:ii$outcomes) {
      got <- genD(function(param) {
        rpf.prob(ii, ii.p, param)[ox]
      }, at, method.args=list(eps=0.01, d=0.01, r=2))
      deriv[,ox] <- got$D
    }
    if (ii$factors == 1) {
      expect_equal(analytic$gradient, c(deriv[1,]), tolerance=1e-6)
      expect_equal(analytic$hessian, c(deriv[2,]), tolerance=1e-3)
    } else {
      expect_equal(analytic$gradient,
                   c(dir %*% deriv[c(1:ii$factors),]), tolerance=1e-6)
      expect_equal(analytic$hessian,
                   c(dir %*% deriv[c(ii$factors+triSize(1:ii$factors)),]), tolerance=1e-3)
    }
  })
}
