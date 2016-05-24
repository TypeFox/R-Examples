
##

context("Taper classes and methods")

taps.o = c(0,1:10,100)
taps = c(0,1:10,100)
ataps <- as.tapers(taps)
ataps.s <- as.tapers(taps, setspan = TRUE)
ms.taps <- minspan(taps)
ms.ataps <- minspan(ataps)

test_that("coercion is functioning",{
  
  expect_is(ataps, 'tapers')
  expect_is(ataps.s, 'tapers')
  
  expect_is(as.vector(ataps), 'integer')
  expect_is(as.vector(ataps.s), 'integer')
  
  expect_is(as.data.frame(ataps), 'data.frame')
  expect_is(as.data.frame(ataps.s), 'data.frame')
  
  expect_is(summary(ataps), 'summary.tapers')
  expect_is(summary(ataps.s), 'summary.tapers')
  
  expect_is(ms.taps, 'integer')
  expect_is(ms.ataps, 'tapers')
  
})

##

context("Taper constraints -- using as.tapers")

test_that("constrained-range is correct",{
  
  expect_equal(min(ataps),  1)
  expect_equal(max(ataps), max(taps))
  
  expect_equal(min(ataps.s), min(minspan(taps)))
  expect_equal(max(ataps.s), max(minspan(taps)))
  
})

##

context("Taper constraints -- constraint algorithms")

test_that("environment variables are protected",{
  
  expect_equal(taps,taps.o)

  expect_is(xx <- ctap_simple(taps), 'integer')
  expect_equal(taps,taps.o)
  
  expect_is(xx <- ctap_simple(ataps), 'tapers')
  expect_equal(taps,taps.o)

  expect_is(xx <- ctap_simple_rcpp(taps), 'integer')
  expect_equal(taps,taps.o)
  
  expect_is(xx <- ctap_simple_rcpp(ataps), 'tapers')
  expect_equal(taps,taps.o)

  expect_warning(xx <- ctap_loess(taps)) # because a sequence is not given
  expect_equal(taps,taps.o)
  
  expect_is(xx <- suppressWarnings(ctap_loess(taps)), 'integer')
  expect_equal(taps,taps.o)
  
  expect_is(xx <- suppressWarnings(ctap_loess(ataps)), 'tapers')
  expect_equal(taps,taps.o)
  
  expect_is(xx <- constrain_tapers(taps, verbose = FALSE), 'integer')
  expect_equal(taps,taps.o)
  
  expect_is(xx <- constrain_tapers(ataps, verbose = FALSE), 'tapers')
  expect_equal(taps,taps.o)
  
})

test_that("constraint coercion is functioning",{
  
  expect_is(constrain_tapers(taps, verbose = FALSE), 'integer')
  expect_is(constrain_tapers(ataps, verbose = FALSE), 'tapers')
  
  expect_is(ctap_simple(taps), 'integer')
  expect_is(ctap_simple(ataps), 'tapers')
  
  expect_is(ctap_simple_rcpp(taps), 'integer')
  expect_is(ctap_simple_rcpp(ataps), 'tapers')
  
  expect_is(suppressWarnings(ctap_loess(taps)), 'integer')
  expect_is(suppressWarnings(ctap_loess(ataps)), 'tapers')
  
})

##

context("Taper constraints -- Rcpp implementaion of ctap_simple")

test_that("constrained-range is correct",{
  
  taps <- c(0,1:10,100)
  
  taps.c <- ctap_simple_rcpp(taps, maxslope=1)
  taps.c2 <- ctap_simple_rcpp(taps, maxslope=2)
  
  expect_equal(min(taps.c), 1)
  expect_equal(min(taps.c2), 1)
  expect_equal(max(taps.c), 11)
  expect_equal(max(taps.c2), 12)
  
})

test_that("bad input is handled correctly",{
  
  expect_equal(ctap_simple_rcpp(NA), 1)
  expect_warning(ctap_simple_rcpp(Inf))
  expect_error(rcpp_ctap_simple(NULL))
  expect_equal(ctap_simple_rcpp(NULL), integer(0))
  expect_error(ctap_simple_rcpp(1, maxslope=-1))
  
})

##

context("Taper constraints -- through minspan")

test_that("Length and positivity requirements are checked correcly",{
  
  expect_error(minspan(1))
  expect_error(minspan(0))
  expect_error(minspan(-1))
  expect_error(minspan(-1:0))
  
})

test_that("the result is limited by section length", {
  
  ms. <- minspan(0:2)
  expect_equal(min(ms.), 1)
  expect_equal(max(ms.), 1)
  
  ms. <- minspan(0:3)
  expect_equal(min(ms.), 1)
  expect_equal(max(ms.), 2)
  
  ms. <- minspan(0:4)
  expect_equal(min(ms.), 1)
  expect_equal(max(ms.), 2)
  
  ms. <- minspan(0:5)
  expect_equal(min(ms.), 1)
  expect_equal(max(ms.), 3)
  
  ms. <- minspan(0:6)
  expect_equal(min(ms.), 1)
  expect_equal(max(ms.), 3)
  
  ms. <- minspan(0:7)
  expect_equal(min(ms.), 1)
  expect_equal(max(ms.), 4)
})

test_that("strange values are dealt with", {
  
  expect_warning(ms. <- minspan(c(0:7,Inf)))
  expect_equal(min(ms.), 1)
  expect_equal(max(ms.), 4)
  expect_equal(length(ms.), 9)
  
  ms. <- minspan(c(0:7,NA))
  expect_equal(min(ms.), 1)
  expect_equal(max(ms.), 4)
  expect_equal(length(ms.), 9)
  
  ms. <- minspan(c(0:7,""))
  expect_equal(min(ms.), 1)
  expect_equal(max(ms.), 4)
  expect_equal(length(ms.), 9)
  
  ms. <- minspan(c(0:7,NULL))
  expect_equal(min(ms.), 1)
  expect_equal(max(ms.), 4)
  expect_equal(length(ms.), 8) # instead of 9
  
})

##
