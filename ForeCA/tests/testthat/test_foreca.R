context("foreca")

set.seed(2)
kNumVariables <- 5
kNumObs <- 200
kSeries <- matrix(arima.sim(n = kNumObs * kNumVariables, list(ar = -0.7)), 
                  ncol = kNumVariables)
kSeries[, 1] <- cumsum(kSeries[, 1])
kSeries <- ts(kSeries)
# rotate series to have non-trivial correlations
kMixing <- matrix(rnorm(kNumVariables^2), ncol = kNumVariables)
kSeries <- ts(kSeries %*% kMixing)

kSeriesCentered <- sweep(kSeries, 2, colMeans(kSeries), "-")

UU <- whiten(kSeries)$U
ww0 <- initialize_weightvector(UU, method = "rnorm")

yy.UU <- UU %*% t(ww0)
yy.Series <- ts(kSeries %*% t(ww0))

context("Foreca: one weightvector")

test_that("only accepts whitened input", {
            expect_error(foreca.one_weightvector(U = kSeries))
          })

one.weight <- foreca.one_weightvector(U = UU)

test_that("Algorithm is not implemented", {
  expect_error(foreca(XX, n.comp = 1, algorithm.type = "foo"))
})


test_that("weightvector has norm 1", {
  expect_equal(base::norm(one.weight$weightvector, "2"), 1)
})

#test_that("Omega from algorithm objective matches Omega estimate of best ForeC", {  
#  omega.direct <- Omega(series = one.weight$score,
#                        spectrum.control = one.weight$spectrum.control,
#                        entropy.control = one.weight$entropy.control)
  
#  expect_equal(one.weight$Omega, omega.direct)
#})

test_that("Omega from best.f spectrum matches reported Omega estimate", {
  omega.spec <- Omega(mvspectrum.output = one.weight$best.f,
                      entropy.control = one.weight$entropy.control)
  expect_equal(one.weight$Omega, c(omega.spec))
})

test_that("score is computed correctly with weightvector", {
  expect_equivalent(ts(UU %*% t(one.weight$weightvector)), one.weight$score)
})

test_that("score has 0 mean and variance 1", {
  expect_equal(0, mean(one.weight$score))
  expect_equal(1, sd(one.weight$score))
})

test_that("spectral density is normalized",{
  expect_equal(0.5, sum(one.weight$best.f))
})


context("Foreca: multiple weightvectors")

test_that("only accepts whitened input", {
  expect_error(foreca.multiple_weightvectors(U = kSeries))
})

kNComp <- 4
mult.weights <- foreca.multiple_weightvectors(U = UU, n.comp = kNComp)
WW <- mult.weights$weightvectors

test_that("weightvectors have right dimension", {
  expect_equal(dim(WW), c(ncol(UU), kNComp))
})

test_that("weightvectors have norm 1", {
  expect_equal(apply(WW, 2, base::norm, "2"), 
               rep(1, kNComp))
})

test_that("weightvectors are orthogonal", {
  expect_equal(diag(1, kNComp), t(WW) %*% WW)
})

test_that("it returns correct number of ForeCs", {
  expect_equal(kNComp, ncol(mult.weights$scores))
})

test_that("ForeCs are whitened", {
  expect_silent(check_whitened(mult.weights$scores, FALSE))
})

names(mult.weights$Omega) <- NULL
test_that("ForeCs are sorted by Omega", {
  expect_equal(kNComp:1, rank(mult.weights$Omega))
})

test_that("ForeC 1 Omega is greater or equal to all input series", {
  expect_true(all(mult.weights$Omega[1] >= mult.weights$Omega.series))
})

test_that("ForeCs scores and series * weightvectors match", {
  expect_equivalent(mult.weights$scores, 
                    ts(UU %*% WW))
})


context("ForeCA: foreca")
mod.foreca <- foreca(kSeries, n.comp = kNComp)

test_that("scores have correct names", {
  expect_equal(colnames(mod.foreca$scores), 
               paste0("ForeC", seq_len(kNComp)))
})

test_that("loadings have correct names", {
  expect_equal(colnames(mod.foreca$loadings), 
               paste0("ForeC", seq_len(kNComp)))
  
  expect_equal(rownames(mod.foreca$loadings), 
               colnames(kSeries))
})

test_that("loadings correctly compute the scores from original series", {
  sl <- ts(scale(kSeries, scale = FALSE) %*% mod.foreca$loadings)
  expect_equal(sl, mod.foreca$scores)
})