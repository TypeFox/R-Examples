
##

context("Parabolic weighting -- Rcpp implementation")

PWS <- function(n, verbose=FALSE){
  n <- as.integer(n)
  pw <- parabolic_weights_rcpp(n)
  pws. <- sum(pw[['taper_weights']])
  if (verbose) message(n," sums to ",pws.)
  return(pws.)
}

test_that("internal checks are working",{
  expect_error(PWS(-1))
  expect_error(PWS(1:2))
})

test_that("weights sum to 1 when n > 0",{
  expect_equal(PWS(0), 0)
  expect_equal(PWS(1), 1)
  expect_equal(PWS(10), 1)
  expect_equal(PWS(100), 1)
  expect_equal(PWS(1e3), 1)
  expect_equal(PWS(1e4), 1)
})

##
