context("mvspectrum")

kNumVariables <- 3
kNumObs <- 1000 - 1
kSeries <- matrix(arima.sim(n = kNumObs * kNumVariables, list(ar = -0.7)),
                  ncol = kNumVariables)
kSeries[, 1] <- cumsum(kSeries[, 1])
kSeriesCentered <- sweep(kSeries, 2, colMeans(kSeries), "-")
whitenedSeries <- whiten(kSeries)$U

num.freqs.exp <- floor(kNumObs /2)

test_that("spec.pgram does not depend on mean", {
  pgram.Series <- spec.pgram(kSeries, plot = FALSE, fast = FALSE)
  pgram.SeriesCentered <- spec.pgram(kSeriesCentered, plot = FALSE, 
                                     fast = FALSE)
  pgram.SeriesCentered$series <- pgram.Series$series
  
  expect_equal(pgram.Series,
               pgram.SeriesCentered)
  
  sapply(names(pgram.Series),
         function(x) {
           expect_equal(pgram.Series[[x]],
                        pgram.SeriesCentered[[x]])
         })
})


for (mm in kMvspectrumMethods) {
  test.msg <- paste0("Testing method ", mm, "\n")
  sc.tmp <- list(method = mm)
  
  spec.Series <- mvspectrum(kSeries, mm)
  spec.SeriesCentered <- mvspectrum(kSeriesCentered, mm)
  
  test_that("mvspectrum is independent of mean", {
    expect_equal(spec.SeriesCentered, spec.Series,
                 info = test.msg)
  })
  
  test_that("mvspectrum has right dimensions", {  
    expect_equal(c(num.freqs.exp, kNumVariables, kNumVariables), 
                 dim(spec.Series),
                 info = test.msg)
  })
  
  test_that("mvspectrum is real valued in diagonal", {
    all.diags <- c(apply(spec.Series, 1, function(x) Im(diag(x))))
    # diagonal is real valued
    expect_equal(rep(0, length(all.diags)),
                 all.diags,
                 info = test.msg)
  })
  
  test_that("mvspectrum is Hermitian for every frequency", {
    freq.conj.diff <- apply(spec.Series, 1, function(x) base::norm(x - Conj(t(x)), "2"))
    # for every frequency is Hermitian
    expect_equal(rep(0, length(freq.conj.diff)),
                 freq.conj.diff,
                 info = test.msg)
  })
  
  test_that("mvspectrum is positive semi-definite for every frequency", {
    # positive semi-definite for each frequency
    lambdas <- apply(spec.Series, 1, function(x) eigen(x)$values)

    lambdas.flat <- c(lambdas)
    expect_equal(Im(lambdas.flat), rep(0, length(lambdas.flat)))
    
    lambdas.flat <- Re(lambdas.flat)
    lambdas.pos <- (round(lambdas.flat, 4) >= 0)
    expect_true(all(lambdas.pos),
                info = paste0(test.msg, ";\n ", 
                              sum(!lambdas.pos), " are negative"))
  })
}


context("Univariate spectra")

for (mm in kMvspectrumMethods) {
  test.msg <- paste0("Testing method ", mm, "\n")
  
  test_that("mvspectrum(x) must give sum close to var(x) / 2", {
    spec.yy.2 <- mvspectrum(scale(kSeries[, 2], scale = FALSE), 
                            method = mm)
    expect_equal(1, var(kSeries[, 2]) / 2 / sum(spec.yy.2), 
                 tolerance = 3 * sd(kSeries[, 2]) / sqrt(nrow(kSeries)),
                 info = test.msg)
  })
    
  if (mm == 'direct') {
    
    test_that("get_spectrum_from_mvspectrum works", {
      test.msg <- paste0("Testing method direct \n")
      
      spec.Series <- mvspectrum(kSeries, "direct")
      spec.Series.1 <- mvspectrum(kSeries[, 1], "direct")
      spec.Series.2 <- mvspectrum(kSeries[, 2], "direct")
      expect_equal(c(spec.Series.1), 
                   c(get_spectrum_from_mvspectrum(spec.Series, 1)),
                   info = test.msg)
      expect_equal(c(spec.Series.2), 
                   c(get_spectrum_from_mvspectrum(spec.Series, 2)),
                   info = test.msg)
      # if which is unspecified return all
      all.spectra <- get_spectrum_from_mvspectrum(spec.Series)
      expect_equal(ncol(spec.Series), ncol(all.spectra),
                   info = test.msg)
    })
    
    test_that("get_spectrum_from_mvspectrum works for univariate spectra", {
      test.msg <- paste0("Testing method direct \n")
      
      spec.Series.1 <- mvspectrum(kSeries[, 1], "direct")
      expect_equal(c(spec.Series.1), 
                   c(get_spectrum_from_mvspectrum(spec.Series.1, 1)),
                   info = test.msg)
      # if which is unspecified return all
      expect_equal(c(spec.Series.1), 
                   c(get_spectrum_from_mvspectrum(spec.Series.1)),
                   info = test.msg)
      expect_error(get_spectrum_from_mvspectrum(spec.Series.1, 2))
    })
    
    beta.tmp <- cbind(rnorm(ncol(kSeries)))
    
    test_that("spectrum_of_linear_combination gives same as direct estimation", {
      yy.tmp <- kSeries %*% beta.tmp
      spec.Series <- mvspectrum(kSeries, method = mm)
      spec.Series.comb <- spectrum_of_linear_combination(spec.Series, 
                                                         beta.tmp)
      spec.yy <- mvspectrum(yy.tmp, method = mm)
      expect_equal(c(spec.yy), c(spec.Series.comb),
                   info = test.msg)
    })
    
  }

  test_that("for basis vector spectrum_of_linear_combination and get_spectrum_from_mspectrum coincide", {
    e.vec <- c(1, rep(0, ncol(kSeries) - 1))
    
    spec.Series <- mvspectrum(whitenedSeries, method = mm, normalize = TRUE)
    spec.Series.comb <- spectrum_of_linear_combination(spec.Series, 
                                                               e.vec)
    spec.Series.get <- get_spectrum_from_mvspectrum(spec.Series, 1)
    expect_equal(spec.Series.comb, spec.Series.get, tolerance = 1e-3,
                 info = test.msg)
  })
}


context("Normalizing spectra\n")

for (mm in kMvspectrumMethods) {
  test.msg <- paste0("Testing method ", mm, "\n")
  
  spec.whitenedSeries <- mvspectrum(whitenedSeries, method = mm)
  spec.Series <- mvspectrum(kSeries, method = mm)
  
  test_that("it throws error if input is not whitened", {    
    expect_error(normalize_mvspectrum(spec.whitenedSeries, Sigma.X = 1),
                 info = test.msg)
    expect_error(normalize_mvspectrum(spec.whitenedSeries, Sigma.X = NA),
                 info = test.msg)
  })
  
  norm.spec.Series <- normalize_mvspectrum(spec.whitenedSeries)
  sum.norm.spec.Series <- apply(norm.spec.Series, 2:3, sum)
  cov.norm.spec.Series <- mvspectrum2wcov(norm.spec.Series)
  off.diag <- sum.norm.spec.Series - diag(diag(sum.norm.spec.Series))

  test_that("mvspectrum has attribute 'normalize' = TRUE only if normalize=TRUE", {
    expect_false(attr(spec.Series, "normalized"),
                 info = test.msg)
    expect_true(attr(norm.spec.Series, "normalized"),
                info = test.msg)
  })
  
  
  test_that("normalize gives 0 real off-diagonals", {
    expect_equal(matrix(0, ncol = ncol(kSeries), nrow = ncol(kSeries)),
                 Re(off.diag),
                 info = test.msg)
  })
  
  test_that("normalize gives Hermitian", {
    expect_equal(off.diag, 
                 Conj(t(off.diag)),
                 info = test.msg)
  })
    
  test_that("normalize has real valued diagonal equal to 0.5", {
    # imaginary part is 0 in diagonal
    expect_equal(rep(0, ncol(kSeries)),
                 Im(diag(sum.norm.spec.Series)),
                 info = test.msg)
    expect_equal(rep(0.5, ncol(kSeries)), 
                 Re(diag(sum.norm.spec.Series)),
                 info = test.msg)
    
    # covariance estimate is identity
    expect_equal(diag(1, ncol(kSeries)), 
                 cov.norm.spec.Series,
                 info = test.msg)
  })
  
  #test_that("normalize gives a matrix that is cov(X) / 2 for non-whitened data", {
  #  norm.spec.Series <- normalize_mvspectrum(spec.Series, 
  #                                           Sigma.X = cov(kSeries))
  #  sum.norm.spec.Series <- apply(norm.spec.Series, 2:3, sum)
  #  sum.norm.spec.Series <- 2 * Re(sum.norm.spec.Series)
  #  expect_true(cor(c(cov(kSeries)), c(sum.norm.spec.Series)) > 0.95,
  #              info = test.msg)
  #})
  
  test_that("normalize makes it add up to 0.5", {
    norm.spec.yy.2 <- mvspectrum(scale(kSeries[, 2]), 
                                 method = mm, normalize = TRUE)
    expect_equal(0.5, sum(norm.spec.yy.2),
                 info = test.msg)
  })
  
  beta.tmp <- t(initialize_weightvector(num.series = ncol(kSeries),
                                        method = 'rnorm'))
  
  spec.whitenedSeries <- mvspectrum(whitenedSeries, method = mm)      
  norm.spec.Series <- normalize_mvspectrum(spec.whitenedSeries)
  
  test_that("if multivariate is normalized, and we use L1 norm vector, so is its linear combination", {
    beta.tmp.norm <- beta.tmp / base::norm(beta.tmp, "2")
    spec.dens.est <- spectrum_of_linear_combination(norm.spec.Series, 
                                                            beta.tmp.norm)
    tmp <- try(check_mvspectrum_normalized(spec.dens.est),
               silent = TRUE)
    expect_false(inherits(tmp, "try-error"))
  })
  
  test_that("if multivariate is normalized, so is each univariate spectra (in the diagonal).", {
    for (ii in seq_len(dim(norm.spec.Series)[2])) {
      spec.dens.est <- get_spectrum_from_mvspectrum(norm.spec.Series, 
                                                    which = ii)
      tmp <- try(check_mvspectrum_normalized(spec.dens.est),
                 silent = TRUE)
      expect_false(inherits(tmp, "try-error"),
                   info = paste0(test.msg, ": series ", ii))
    }
  })
  
  # Not true anymore if normalization is done by pre-multiplying
  # by inverse of sum.
  test_that("normalized_mvspectrum is real valued in diagonal", {
    all.diags <- c(apply(norm.spec.Series, 1, function(x) Im(diag(x))))
    # diagonal is real valued
    expect_equal(rep(0, length(all.diags)),
                 all.diags,
                 info = test.msg)
  })
  
  test_that("normalized_mvspectrum is Hermitian for every frequency", {
    freq.conj.diff <- apply(norm.spec.Series, 1, 
                            function(x) base::norm(x - Conj(t(x)), "2"))
    # for every frequency is Hermitian
    expect_equal(rep(0, length(freq.conj.diff)),
                 freq.conj.diff,
                 info = test.msg)
  })
  

  test_that("normalize mvspectrum is positive semi-definite for every frequency", {
    # positive semi-definite for each frequency
    lambdas <- apply(norm.spec.Series, 1, function(x) eigen(x)$values)
    
    lambdas.flat <- c(lambdas)
    expect_equal(Im(lambdas.flat), rep(0, length(lambdas.flat)))
    
    lambdas.flat <- Re(lambdas.flat)
    lambdas.pos <- (round(lambdas.flat, 4) >= 0)
    expect_true(all(lambdas.pos),
                info = paste0(test.msg, ";\n ", 
                              sum(!lambdas.pos), " are negative"))
  })
    
  test_that("L2 norm = 1 combination of whitened series: spectrum_of_linear_combination", {
    yy.tmp <- whitenedSeries %*% beta.tmp
    
    spec.Series <- mvspectrum(whitenedSeries, method = mm, normalize = TRUE)
    spec.Series.comb <- spectrum_of_linear_combination(spec.Series, 
                                                       beta.tmp)
    expect_equal(0.5, sum(spec.Series.comb))
    spec.yy <- mvspectrum(yy.tmp, method = mm, normalize = TRUE)
    
    #layout(matrix(1:4, ncol = 2))
    #matplot(cbind(spec.Series.comb, spec.yy))
    #plot(spec.Series.comb, spec.yy)
    #plot(log(spec.Series.comb + median(spec.yy)), log(spec.yy + median(spec.yy)))
    expect_true(cor(spec.Series.comb, spec.yy) > 0.95,
                info = test.msg)
  })

  
}
