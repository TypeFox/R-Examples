library(DNAprofiles)
context("Sampling profiles")

test_that(desc = "Sampling profiles",{    
  data(freqsNLngm,freqsNLsgmplus)
  
  # see if sampling profiles works
  n <- 1e3
  x <- sample.profiles(N = n,freqs = freqsNLngm)
  expect_equal(nrow(x),n)
  expect_equal(ncol(x)/2,length(freqsNLngm))
  expect_false(any(is.na(x)))
  
  # with theta>0
  x <- sample.profiles(N = n,freqs = freqsNLngm,theta = 0.1)
  expect_equal(nrow(x),n)
  expect_equal(ncol(x)/2,length(freqsNLngm))
  expect_false(any(is.na(x)))
  
  # with theta=1
  x <- sample.profiles(N = n,freqs = freqsNLngm,theta = 1)
  expect_equal(nrow(x),n)
  expect_equal(ncol(x)/2,length(freqsNLngm))
  expect_false(any(is.na(x)))
  expect_true(all(homozygous(x)))  
  
  # see if sampling profiles works with subset of the markers
  x <- sample.profiles(N = n,markers = names(freqsNLsgmplus),freqs = freqsNLngm)
  expect_equal(nrow(x),n)
  expect_equal(ncol(x)/2,length(freqsNLsgmplus))
  expect_false(any(is.na(x)))
  
  # with theta>0 with subset of the markers
  x <- sample.profiles(N = n, names(freqsNLsgmplus), freqs = freqsNLngm,theta = 0.1)
  expect_equal(nrow(x),n)
  expect_equal(ncol(x)/2,length(freqsNLsgmplus))
  expect_false(any(is.na(x)))
  
  # with theta=1 with subset of the markers
  x <- sample.profiles(N = n, names(freqsNLsgmplus), freqs = freqsNLngm,theta = 1)
  expect_equal(nrow(x),n)
  expect_equal(ncol(x)/2,length(freqsNLsgmplus))
  expect_false(any(is.na(x)))
  expect_true(all(homozygous(x)))  
})



