context("\ninitialize_weightvector")

kNumVariables <- 4
kAvailableMethods <- c("SFA", "SFA.slow", "SFA.fast", 
                       "PCA", "PCA.small", "PCA.large",
                       "rcauchy", "runif", "rnorm", "max")

kNumObs <- 100
kTimeSeries <- matrix(rnorm(kNumVariables * kNumObs), ncol = kNumVariables)
# this makes the first variable the most forecastable one and also the slowes
kTimeSeries[, 1] <- arima.sim(kNumObs, model = list(ar = 0.9))
# is faster than white noise, yet less forecastable than the first series
kTimeSeries[, 2] <- arima.sim(kNumObs, model = list(ar = -0.5))

kWhitened <- whiten(kTimeSeries)$U
f.U.tmp <- mvspectrum(kWhitened, "wosa", normalize = TRUE)

for (mm in kAvailableMethods) {
  test.msg <- paste0("Testing method ", mm, "\n")
  
  test_that("initialize_weightvector has correct length and norm 1", {
    ww.tmp <- initialize_weightvector(kWhitened, f.U.tmp, method = mm)
    expect_equal(kNumVariables, length(ww.tmp))
    # L2 norm 1
    expect_equal(base::norm(ww.tmp, "2"), 1,
                 info = test.msg)
  })
}

test_that("Max is actually max Omega", {
  ww.tmp <- initialize_weightvector(U = kWhitened, f.U = f.U.tmp, method = "max")
  expect_equal(which.max(Omega(kWhitened, spectrum.control = list(method = "wosa"))), 
               which.max(ww.tmp),
               info = test.msg)
  # the first one is the AR(1)
  expect_equal(which.max(Omega(kWhitened, spectrum.control = list(method = "wosa"))), 1,
               info = test.msg)
})

test_that("SFA gives slowest/fastest AR(1) signal: initialize_weightvector", {
  ww.tmp <- initialize_weightvector(U = kWhitened, f.U = f.U.tmp, method = "SFA.slow")
  # the first signal is the positive lag 1 autocorrelation AR(1)
  expect_equal(which.max(abs(ww.tmp)), 1,
               info = test.msg)
  ww.tmp <- initialize_weightvector(U = kWhitened, f.U = f.U.tmp, method = "SFA.fast")
  # the fastest signal is the negative lag 1 autocorrelation AR(1)
  expect_equal(which.max(abs(ww.tmp)), 2,
               info = test.msg)
  ww.tmp <- initialize_weightvector(U = kWhitened, f.U = f.U.tmp, method = "SFA")
  # the first signal is more forecastable than the 2nd one (0.9 vs -0.5 lag 1 autocorrelation)
  expect_equal(which.max(abs(ww.tmp)), 1,
               info = test.msg)
})
