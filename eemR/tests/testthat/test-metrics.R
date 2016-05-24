context("metrics")

file <- system.file("extdata/cary/scans_day_1/sample1.csv", package = "eemR")
eems <- eem_read(file)

test_that("Test Cobble's peaks", {

  metrics <- eem_coble_peaks(eems, verbose = FALSE)

  b <- 1.5452981
  t <- 1.060331225
  a <- 3.731835842 # approximative value since some wl do not exist
  m <- 2.443755269 # approximative value since some wl do not exist
  c <- 1.815422177 # approximative value since some wl do not exist

  expect_equal(b, metrics$b)
  expect_equal(t, metrics$t)
  expect_equal(a, metrics$a)
  expect_equal(m, metrics$m, tolerance = 0.01)
  expect_equal(c, metrics$c, tolerance = 0.001)

})

test_that("Test fluorescence index (FI)", {

  metrics <- eem_fluorescence_index(eems, verbose = FALSE)

  fi <- 1.124091625 / 0.888762951

  expect_equal(fi, metrics$fi)

})

test_that("Test biological fluorescence index (BIX)", {

  metrics <- eem_biological_index(eems, verbose = FALSE)

  bix <- 1.735240459 / 2.456928968

  expect_equal(bix, metrics$bix)

})




