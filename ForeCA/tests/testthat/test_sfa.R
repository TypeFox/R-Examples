context("foreca")

set.seed(2)
kNumVariables <- 4
kNumObs <- 200
kSeries <- matrix(arima.sim(n = kNumObs * kNumVariables, list(ar = -0.6)), 
                  ncol = kNumVariables)
kSeries[, 1] <- cumsum(kSeries[, 1])
# rotate series to have non-trivial correlations
kMixing <- matrix(rnorm(kNumVariables^2), ncol = kNumVariables)
kSeries <- kSeries %*% kMixing

sfa.est <- sfa(kSeries)

rSFA.installed <- requireNamespace("rSFA", quietly = TRUE)

if (rSFA.installed) {
  rsfa.est <- rSFA::sfa1(kSeries)
}

.aux_eta <- function(x) {
  return(mean(diff(scale(x))^2))
}

# remove names 
attr(sfa.est$scores, "dimnames") <- NULL

test_that("sfa() returns signals in increasing order", {
  etas <- apply(sfa.est$scores, 2, .aux_eta)
  rhos <- apply(sfa.est$scores, 2, function(x) acf(x, plot = FALSE)$acf[2])
  
  expect_equal(order(rhos, decreasing = TRUE), seq_len(kNumVariables))
  expect_equal(order(etas), seq_len(kNumVariables))
})

test_that("sfa() returns orthonormal signals with mean 0", {
  expect_equal(cov(sfa.est$scores), 
               diag(1, kNumVariables),
               tol = 1e-5)
  expect_equal(colMeans(sfa.est$scores),
               rep(0, kNumVariables),
               tol = 1e-5)
})


if (rSFA.installed) {

  test_that("rsfa() returns signals in increasing order", {
    etas <- apply(rsfa.est$y, 2, .aux_eta)
    rhos <- apply(rsfa.est$y, 2, function(x) acf(x, plot = FALSE)$acf[2])
    
    expect_equal(order(rhos, decreasing = TRUE), seq_len(kNumVariables))
    expect_equal(order(etas), seq_len(kNumVariables))
  })
  
  test_that("sfa() gives the same as sfa1() from rSFA package", {
     expect_equal(abs(sfa.est$scores), abs(rsfa.est$y), tol = 1e-2,
                  info = "extracted slow features are not the same")
     
     loadings.tmp <- sfa.est$loadings
     class(loadings.tmp) <- "matrix"
     dimnames(loadings.tmp) <- NULL
     
     expect_equal(abs(loadings.tmp), abs(t(rsfa.est$SF)), tol = 1e-3)
     expect_equal(c(rsfa.est$avg0), c(sfa.est$center))
  })
}

