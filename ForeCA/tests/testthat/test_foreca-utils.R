context("foreca-utils (S3 methods)")

set.seed(2)
kNumVariables <- 5
kNumObs <- 200
kSeries <- matrix(arima.sim(n = kNumObs * kNumVariables, list(ar = -0.7)), 
                  ncol = kNumVariables)
kSeries[, 1] <- cumsum(kSeries[, 1])
# rotate series to have non-trivial correlations
kMixing <- matrix(rnorm(kNumVariables^2), ncol = kNumVariables)
kSeries <- kSeries %*% kMixing

kSeriesCentered <- sweep(kSeries, 2, colMeans(kSeries), "-")

UU <- whiten(kSeries)$U
ww0 <- initialize_weightvector(UU, method = "rnorm")

yy.UU <- UU %*% t(ww0)
yy.Series <- kSeries %*% t(ww0)

context("S3 methods for one weightvector")

one.weight <- foreca.one_weightvector(U = UU)

test_that("print method works", {
  tmp <- try(print(one.weight),
             silent = TRUE)
  expect_false(inherits(tmp, "try-error"))
})

test_that("summary method works", {
  tmp <- try(summary(one.weight),
             silent = TRUE)
  expect_false(inherits(tmp, "try-error"))
})

test_that("summary method works", {
  tmp <- try(plot(one.weight),
             silent = TRUE)
  expect_false(inherits(tmp, "try-error"))
})


context("S3 methods for foreca")

mod.foreca <- foreca(series = kSeries, n.comp = 4)

test_that("print method works", {
  tmp <- try(print(mod.foreca),
             silent = TRUE)
  expect_false(inherits(tmp, "try-error"))
})

test_that("summary method works", {
  tmp <- try(summary(mod.foreca),
             silent = TRUE)
  expect_false(inherits(tmp, "try-error"))
})

test_that("summary method works", {
  tmp <- try(plot(mod.foreca),
             silent = TRUE)
  expect_false(inherits(tmp, "try-error"))
})


