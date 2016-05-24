context("Omega")

kNumVariables <- 5
kNumObs <- 1000

kTimeSeries <- matrix(rnorm(kNumVariables * kNumObs), ncol = kNumVariables)
# this makes the first variable the most forecastable one and also the slowes
kTimeSeries[, 1] <- arima.sim(kNumObs, model = list(ar = 0.9))
# is faster than white noisel, yet less forecastable than the first series
kTimeSeries[, 2] <- arima.sim(kNumObs, model = list(ar = -0.5))
omegas <- Omega(kTimeSeries)

for (mm in kMvspectrumMethods) {
  test.msg <- paste0("Test method ", mm, "\n")
  sc.tmp <- list(method = mm)

  test_that("returns the right number of values and in the right range", {
    expect_equal(length(omegas), ncol(kTimeSeries),
                info = test.msg)
    expect_true(omegas > 0 && omegas < 100,
                info = test.msg)
  }) 
  
  test_that("AR(1) with positive autocorrelation has the highest foreastability", {
    # AR(1) with 0.9 has highest forecastability
    expect_equal(which.max(omegas), 1)
  })

  if (mm != "direct") {
    test_that("Non-parametric estimates have at most forecastability of direct", {
      omegas.direct <- Omega(kTimeSeries, spectrum.control = list(method = "direct"))
      omegas <- Omega(kTimeSeries, spectrum.control = sc.tmp)
      expect_true(all(round(omegas, 4) <= round(omegas.direct, 4)),
                  info = test.msg)
    })
  }

  test_that("Using prior uniform smoothing gives less forecastability", {
    omegas.wo.prior <- Omega(kTimeSeries[, 1], spectrum.control = sc.tmp)
    omegas.w.prior <- Omega(kTimeSeries[, 1], 
                            spectrum.control = sc.tmp,
                            entropy.control = list(prior.weight = 0.1))
    expect_true(all(omegas.w.prior < omegas.wo.prior),
                info = test.msg)
  })
  
  test_that("White noise is not forecastable", {
    omega.wn <- Omega(rnorm(1e4), spectrum.control = sc.tmp,
                      entropy.control = list(base = 2))
    expect_true(omega.wn < 10,
                info = test.msg)
  })
  
  
  spec.Series.1 <- mvspectrum(kSeries[, 1], method = sc.tmp$method)
  omega.direct <- Omega(kSeries[, 1], spectrum.control = sc.tmp)
  omega.spectrum.estimate <- Omega(mvspectrum.output = spec.Series.1)
  
  test_that("using series or spectrum estimate directly does not make a difference", {
    expect_equal(omega.direct, omega.spectrum.estimate,
                 info = test.msg)
  })
  
  spec.ent.direct <- spectral_entropy(kSeries[, 1], spectrum.control = sc.tmp)
  
  test_that("Omega equals (1 - spectral entropy) * 100", {
    expect_equal(as.numeric(omega.direct), 
                 as.numeric((1 - spec.ent.direct) * 100),
                 info = test.msg)
  })
  
  
}
