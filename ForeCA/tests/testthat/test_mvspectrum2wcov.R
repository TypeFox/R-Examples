

context("\nMvspectrum2wcov computation")

for (mm in kMvspectrumMethods) {
  test.msg <- paste0("Test method ", mm, "\n")
  spec.Series <- mvspectrum(kSeries, mm)
  sigma.hat.freq <- mvspectrum2wcov(spec.Series)
  
  test_that("mvspectrum2wcov does actually compute the cov(X)", {
    # check difference of matrices (use inverse covariance matrix to put them on comparable scale)
    expect_true(base::norm(solve(Sigma, sigma.hat.freq - Sigma), "2") < 0.05,
                info = test.msg)
  })
  
  test_that("it is symmetric for any weights (positive or negative); positive definite for pos weights", {
    sigma.hat.pos.weights <- mvspectrum2wcov(spec.Series, 
                                             kernel.weights = rexp(dim(spec.Series)[1]))
    sigma.hat.general.weights <- mvspectrum2wcov(spec.Series, 
                                                 kernel.weights = rnorm(dim(spec.Series)[1]))
    # symmetric
    expect_equal(sigma.hat.general.weights, t(sigma.hat.general.weights),
                 info = test.msg)
    expect_equal(sigma.hat.pos.weights, t(sigma.hat.pos.weights),
                 info = test.msg)
    # positive semidefinite for positive weights
    expect_true(all(eigen(sigma.hat.pos.weights, only.values = TRUE)$values >= 0),
                info = test.msg)
  })
}


context("\n weightvector2entropy_wcov computation")

for (mm in kMvspectrumMethods) {
  test.msg <- paste0("Test method ", mm, "\n")
  
  spec.Series <- mvspectrum(kSeries, mm)
  sigma.hat.freq <- mvspectrum2wcov(spec.Series)
  
  norm.spec.Series <- mvspectrum(kWhitenedSeries, mm, normalize = TRUE)
  sigma.hat.freq <- mvspectrum2wcov(spec.Series)
  
  yy.norm.spec <- spectrum_of_linear_combination(norm.spec.Series,
                                                         ww.tmp)
  
  int.mvspec <- weightvector2entropy_wcov(ww.tmp,
                                          norm.spec.Series)
  int.mvspec.direct <- weightvector2entropy_wcov(NULL,
                                                 norm.spec.Series,
                                                 yy.norm.spec)
  
  test_that("weightvector and direct spectrum.estimate give the same entropy", {
    expect_equal(int.mvspec, int.mvspec.direct,
                 info = test.msg, tol = 1e-6)
  })
  
  yy.entropy.by.int <- quadratic_form(int.mvspec, ww.tmp)
  yy.entropy.direct <- spectral_entropy(mvspectrum.output = yy.norm.spec)
  
  test_that("spectral entropy by direct estimation is the same as by quadratic form on integrated mvspectrum", {
    expect_equal(yy.entropy.by.int, 
                 as.numeric(yy.entropy.direct),
                 tol = 1e-3,
                 info = test.msg)
  })


  int.mvspec.prior <- 
    weightvector2entropy_wcov(ww.tmp,
                              norm.spec.Series,
                              entropy.control = ec.tmp)
  yy.entropy.prior <- spectral_entropy(mvspectrum.output = yy.norm.spec,
                                       entropy.control = ec.tmp)
  
  test_that("prior spectral entropy is larger than normal spectral entropy", {
    expect_true(yy.entropy.prior > yy.entropy.direct,
                info = test.msg)
  })
  
  yy.entropy.by.int.prior <- quadratic_form(int.mvspec.prior, ww.tmp)
  
  # test_that("spectral entropy by direct estimation is the same as by quadratic form on integrated mvspectrum", {
  #  expect_equal(yy.entropy.by.int.prior,
  #               as.numeric(yy.entropy.prior),
  #               info = test.msg)
  #})
}
