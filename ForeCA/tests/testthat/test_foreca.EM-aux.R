context("foreca.EM-aux")
kNumVariables <- 3
kNumObs <- 500
kSeries <- matrix(arima.sim(n = kNumObs * kNumVariables, list(ar = 0.5)) + 10, 
                  ncol = kNumVariables)
kSeries[, 1] <- cumsum(kSeries[, 1])
kSeriesCentered <- sweep(kSeries, 2, colMeans(kSeries), "-")

UU <- whiten(kSeries)$U

ww0 <- initialize_weightvector(num.series = ncol(UU), method = 'rnorm')
yy.UU <- UU %*% t(ww0)
attr(yy.UU, "whitened") <- TRUE

yy.Series <- kSeries %*% t(ww0)

kSpectrumMethods <- c("direct", "wosa", "multitaper", "mvspec", "pgram")

context("foreca.EM.E_step")
for (mm in kSpectrumMethods) {
  test.msg <- paste0("Test method ", mm, "\n")
  
  test_that("E-step throws error if spectrum is not normalized",{
    # must be normalized spectrum
    expect_error(foreca.EM.E_step(mvspectrum(kSeries, method = mm, 
                                             normalize = FALSE), ww0),
                 info = test.msg)
  })
  
  spec.UU <- mvspectrum(UU, method = mm, normalize = TRUE)
  spec.UU.e_step <- foreca.EM.E_step(spec.UU, ww0)
  spec.yy.UU <- mvspectrum(yy.UU, method = mm, normalize = TRUE)
  
  test_that("E-step returns normalized spectrum", {
    expect_true(all(spec.UU.e_step >= 0))
    expect_equal(sum(spec.UU.e_step), expected = 0.5,
                 info = test.msg)
  })
  
  test_that("E-step computes correct linear combination spectrum", {
    # check that combination has same spectrum as yy
    # same for center vs non-centered data
    avg.spec <- mean(spec.UU.e_step)
    expect_gt(cor(log(spec.UU.e_step + avg.spec), 
                      log(spec.yy.UU + avg.spec)),
              0.8)#, info = test.msg)
  })
}

context("foreca.EM.M_step")
for (mm in kSpectrumMethods) {
  test.msg <- paste0("Test method ", mm, "\n")
  
  test_that("M-step maximizes h", {
    spec.UU <- mvspectrum(UU, mm, normalize = TRUE)
    spec.yy.0 <- foreca.EM.E_step(spec.UU, ww0)
    min.val <- foreca.EM.M_step(spec.UU, spec.yy.0, minimize = TRUE)
    max.val <- foreca.EM.M_step(spec.UU, spec.yy.0, minimize = FALSE)
    # L2 norm = 1
    expect_equal(1, base::norm(min.val$vector, "2"),
                 info = test.msg)
    expect_equal(1, base::norm(max.val$vector, "2"),
                 info = test.msg)
    # pos eigenvalue
    expect_true(min.val$value > 0,
                info = test.msg)
    expect_true(max.val$value > 0,
                info = test.msg)
    
    ww1 <- min.val$vector
    yy.1 <- UU %*% ww1
    spec.yy.1 <- foreca.EM.E_step(spec.UU, ww1)
    expect_true(Omega(mvspectrum.output = spec.yy.0) <= 
                  Omega(mvspectrum.output = spec.yy.1),
                info = test.msg)
  })
}


context("foreca.EM.E_and_M_step")
for (mm in kSpectrumMethods) {
  test.msg <- paste0("Test method ", mm, "\n")
  spec.UU <- mvspectrum(UU, mm, normalize = TRUE)
  spec.yy.0 <- foreca.EM.E_step(spec.UU, ww0)
  
  test_that("E_and_M_step wrapper is the same as E and M step separately", {
    spec.UU.e.step <- foreca.EM.E_step(spec.UU, ww0)
    spec.UU.m.step <- foreca.EM.M_step(spec.UU, spec.UU.e.step)

    expect_equal(spec.UU.m.step$vector, 
                 foreca.EM.E_and_M_step(ww0, spec.UU),
                 info = test.msg)
  })
}

context("foreca.EM.h")
for (mm in kSpectrumMethods) {
  test.msg <- paste0("Test method ", mm, "\n")
  spec.UU <- mvspectrum(UU, mm, normalize = TRUE)
  spec.yy.0 <- foreca.EM.E_step(spec.UU, ww0)
  
  test_that("h and Omega add up to 1", {
    expect_equal(object = Omega(mvspectrum.output = spec.yy.0) / 100 + 
                    foreca.EM.h(ww0, spec.UU),
                 expected = 1,
                 check.names = FALSE,
                 check.attributes = FALSE,
                 tol = 1e-3,
                 info = test.msg)
  })

  spec.ent.yy.0 <- foreca.EM.h(ww0, spec.UU)
  one.step <- foreca.EM.M_step(spec.UU, spec.yy.0, 
                               entropy.control = list(prior.weight = 0))
  ww1 <- one.step$vector
  spec.yy.1 <- foreca.EM.E_step(spec.UU, ww1)
  spec.ent.yy.1 <- foreca.EM.h(ww1, spec.UU)

  # bound on spectral entropy
  spec.ent.bound.yy.1 <- foreca.EM.h(ww1, spec.UU, ww0)
  
  test_that("eigen value is close to bound", {
    expect_equal(object = one.step$value,
                 expected = spec.ent.bound.yy.1,
                 tol = 1e-2,
                 info = test.msg)
  })  
  
  test_that("M-step minimizes entropy bound via eigen value inequality", {
    expect_lt(object = spec.ent.bound.yy.1, # min eigenvalue inequality
              expected = spec.ent.yy.0)#,  # iteration 0
              # info = test.msg)
  })
  
  test_that("One more minimization via KL divergence inequality", {
    expect_lt(object = spec.ent.yy.1,# KL divergence inequality
              expected = spec.ent.bound.yy.1)#,  # iteration 0
              # info = test.msg)
  })  
}

