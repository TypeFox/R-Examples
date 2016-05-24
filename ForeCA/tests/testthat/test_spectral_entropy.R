
kNumVariables <- 5
kNumObs <- 1000

kTimeSeries <- matrix(rnorm(kNumVariables * kNumObs), ncol = kNumVariables)
# this makes the first variable the most forecastable one and also the slowes
kTimeSeries[, 1] <- arima.sim(kNumObs, model = list(ar = 0.9))
# is faster than white noisel, yet less forecastable than the first series
kTimeSeries[, 2] <- arima.sim(kNumObs, model = list(ar = -0.5))

context("spectral_entropy")

test_that("spectral_entropy is equal to 1 for perfect uniform spectrum", {
  
  unif.spec.ent <- spectral_entropy(mvspectrum.output = rep(0.1, 10))
  expect_equal(1, c(unif.spec.ent))
})

for (mm in kMvspectrumMethods) {
  test.msg <- paste0("Test method ", mm, "\n")
  sc.tmp <- list(method = mm)
  
  test_that("spectral_entropy is between 0 and 1", {
    SE <- apply(kTimeSeries, 2, spectral_entropy, 
                spectrum.control = sc.tmp)
    
    expect_true(all(SE > 0 && SE < 1),
                info = test.msg)
    # AR(1) with 0.9 has lowest entropy forecastability
    expect_equal(which.min(SE), 1,
                 info = test.msg)
  })
  

  
  test_that("Independent of location/scale: spectral_entropy", {
    se.orig <- spectral_entropy(kTimeSeries[, 1], spectrum.control = sc.tmp)
    se.center <- spectral_entropy(kTimeSeries[, 1] - mean(kTimeSeries[, 1]), 
                                  spectrum.control = sc.tmp)
    se.center.scale <- spectral_entropy(scale(kTimeSeries[, 1]), 
                                        spectrum.control = sc.tmp)
    expect_equal(se.orig, se.center,
                 info = test.msg)
    expect_equal(se.center, se.center.scale,
                 info = test.msg)
  })
  
  test_that("White noise is has high spectral entropy", {
    se.wn <- spectral_entropy(rnorm(1e4), spectrum.control = sc.tmp)
    expect_true(se.wn > 0.9, info = test.msg)
  })

  test_that("Prior uniform increases spectral entropy", {
    se.wn <- spectral_entropy(rnorm(1e4), spectrum.control = sc.tmp)
    expect_true(se.wn > 0.9,
                info = test.msg)
  })

  xx.tmp <- kTimeSeries[, 1]
  eps <- rnorm(length(xx.tmp))
  test_that("signal + noise has larger entropy than signal", {
    expect_true(spectral_entropy(xx.tmp, 
                                 spectrum.control = sc.tmp) <
                spectral_entropy(xx.tmp + eps, 
                                 spectrum.control = sc.tmp),
                info = test.msg)
  })
  
  test_that("Adding prior weight increses spectral entropy", {
    se.wo.prior <- spectral_entropy(xx.tmp, 
                                    spectrum.control = sc.tmp)
    se.w.prior <- spectral_entropy(xx.tmp, 
                                   spectrum.control = sc.tmp,
                                   entropy.control = list(prior.weight = 0.1))
    expect_true(se.w.prior > se.wo.prior,
                info = test.msg)
  })
  
  xx.spec <- mvspectrum(xx.tmp, method = sc.tmp$method)
  spec.ent.direct <- spectral_entropy(xx.tmp, spectrum.control = sc.tmp)
  spec.ent.spectrum.estimate <- spectral_entropy(mvspectrum.output = xx.spec,
                                                 spectrum.control = sc.tmp)
  
  test_that("estimating directly is the same as providing the spectrum.estimate", {
     expect_equal(spec.ent.direct, spec.ent.spectrum.estimate,
                  info = test.msg)
  })
}