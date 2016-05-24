
##

context("Spectrum estimation tools")

tol <- 0.07

test_that("classes are correct",{
  
  set.seed(1234)
  x <- rnorm(100)
  pd <- spectrum(x, plot=FALSE)
  pc <- psdcore(x, plot = FALSE, verbose = FALSE)
  pa <- pspectrum(x, plot = FALSE, verbose = FALSE)
  
  expect_is(pd, 'spec')
  expect_is(pc, c('amt','spec'))
  expect_is(pa, c('amt','spec'))
  
  expect_is(normalize(pd, verbose = FALSE), 'spec')
  expect_is(normalize(pa, verbose = FALSE), c('amt','spec'))
  
  expect_is(pd[['taper']], 'numeric')
  expect_is(pa[['taper']], 'tapers')
  
})

test_that("pspectrum results are accurate",{
  
  set.seed(1234)
  x <- rnorm(100)
  varx <- var(x)
  twovar <- 2*varx
  
  xt <- ts(x, frequency=1)
  xt2 <- ts(x, frequency=10)
  
  pc <- pspectrum(xt, plot = FALSE, verbose = FALSE)
  pc2 <- pspectrum(xt2, plot = FALSE, verbose = FALSE)
  
  # make sure Nyquist frequencies are correct
  fn <- max(pc[['freq']])
  fn2 <- max(pc2[['freq']])
  expect_equal(fn, frequency(xt)/2)
  expect_equal(fn2, frequency(xt2)/2)
  
  # and normalization is for a single-sided spectrum
  nf <- pc[['numfreq']]
  nf2 <- pc2[['numfreq']]
  psum <- sum(pc[['spec']])
  psum2 <- sum(pc2[['spec']])
  expect_equal(fn*psum/nf, varx, tolerance=tol)
  expect_equal(fn2*psum2/nf2, varx, tolerance=tol)
  
  # normalization effects
  pcn <- normalize(pc, verbose = FALSE)
  pcn2 <- normalize(pc2, verbose = FALSE)
  nnf <- pcn[['numfreq']]
  nnf2 <- pcn2[['numfreq']]
  psumn <- sum(pcn[['spec']])
  psumn2 <- sum(pcn2[['spec']])
  expect_equal(fn*psumn/nnf, varx, tolerance=tol)
  expect_equal(fn2*psumn2/nnf2, varx, tolerance=tol)
  
})

test_that("psdcore results are accurate",{
  
  set.seed(1234)
  x <- rnorm(100)
  varx <- var(x)
  twovar <- 2*varx
  
  xt <- ts(x, frequency=1)
  xt2 <- ts(x, frequency=10)
  
  pc <- psdcore(xt, plot = FALSE, verbose = FALSE)
  pc2 <- psdcore(xt2, plot = FALSE, verbose = FALSE)
  
  # make sure Nyquist frequencies are correct
  fn <- max(pc[['freq']])
  fn2 <- max(pc2[['freq']])
  expect_equal(fn, frequency(xt)/2)
  expect_equal(fn2, frequency(xt2)/2)
  
  # and normalization is for a single-sided spectrum
  nf <- pc[['numfreq']]
  nf2 <- pc2[['numfreq']]
  psum <- sum(pc[['spec']])
  psum2 <- sum(pc2[['spec']])
  expect_equal(psum/nf, twovar, tolerance=tol)
  expect_equal(psum2/nf2, twovar, tolerance=tol)
  
})

test_that("prewhiten results are accurate",{
  
  set.seed(1234)
  x <- rnorm(100)
  varx <- var(x)
  twovar <- 2*varx
  
  xt <- ts(x, frequency=1)
  xt2 <- ts(x, frequency=10)
  
  pc <- prewhiten(xt, plot = FALSE, verbose = FALSE)
  pc2 <- prewhiten(xt2, plot = FALSE, verbose = FALSE)
  
  # make sure Nyquist frequencies are correct
  fn <- frequency(pc[['prew_lm']])
  fn2 <- frequency(pc2[['prew_lm']])
  expect_equal(fn, frequency(xt))
  expect_equal(fn2, frequency(xt2))
  
  pa <- prewhiten(xt, AR.max = 10, plot = FALSE, verbose = FALSE)
  pa2 <- prewhiten(xt2, AR.max = 10, plot = FALSE, verbose = FALSE)
  
  # make sure Nyquist frequencies are correct
  fn <- frequency(pa[['prew_ar']])
  fn2 <- frequency(pa2[['prew_ar']])
  expect_equal(fn, frequency(xt))
  expect_equal(fn2, frequency(xt2))
  
})

test_that("pilot_spec results are accurate",{
  
  set.seed(1234)
  x <- rnorm(100)
  varx <- var(x)
  twovar <- 2*varx
  
  xt <- ts(x, frequency=1)
  xt2 <- ts(x, frequency=10)
  
  pc <- pilot_spec(xt, plot = FALSE, verbose = FALSE)
  pc2 <- pilot_spec(xt2, plot = FALSE, verbose = FALSE)
  
  # make sure Nyquist frequencies are correct
  fn <- max(pc[['freq']])
  fn2 <- max(pc2[['freq']])

  expect_equal(fn, frequency(xt)/2)
  expect_equal(fn2, frequency(xt2)/2)
  
  # and normalization is for a single-sided spectrum
  nf <- pc[['numfreq']]
  nf2 <- pc2[['numfreq']]
  psum <- sum(pc[['spec']])
  psum2 <- sum(pc2[['spec']])
  expect_equal(psum/nf, twovar, tolerance=tol)
  expect_equal(psum2/nf2, twovar, tolerance=tol)
  
  expect_warning(pa <- pilot_spec(xt, remove.AR = TRUE, plot = FALSE, verbose = FALSE)) # because there is no AR structure!
  expect_warning(pa2 <- pilot_spec(xt2, remove.AR = TRUE, plot = FALSE, verbose = FALSE))
  
  # make sure Nyquist frequencies are correct
  fn <- max(pa[['freq']])
  fn2 <- max(pa2[['freq']])
  
  expect_equal(fn, frequency(xt)/2)
  expect_equal(fn2, frequency(xt2)/2)
  
})

##
