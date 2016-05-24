library(DNAprofiles)
context("Sampling relatives")

test_that(desc = "Sampling relatives",{    
  data(freqsNLngm,freqsNLsgmplus)
  
  # sampling parent offspring of single x
  n <- 1e3
  x0 <- sample.profiles(N = 1,freqs = freqsNLngm)
  x <- sample.relatives(x = x0,N = n,type = "PO")
  expect_equal(nrow(x),n)
  expect_equal(ncol(x)/2,length(freqsNLngm))
  expect_false(any(is.na(x)))
  
  # with theta>0
  x0 <- sample.profiles(N = 1,freqs = freqsNLngm)
  x <- sample.relatives(x = x0,N = n,type = "PO",theta = 0.1)
  expect_equal(nrow(x),n)
  expect_equal(ncol(x)/2,length(freqsNLngm))
  expect_false(any(is.na(x)))
  
  # with theta=1
  x0 <- sample.profiles(N = 1,freqs = freqsNLngm)
  x <- sample.relatives(x = x0,N = n,type = "PO",theta = 1)  
  expect_equal(nrow(x),n)
  expect_equal(ncol(x)/2,length(freqsNLngm))
  expect_false(any(is.na(x)))   
  
  # same for full sibs with subset of markers
  x <- sample.relatives(x = x0,N = n,markers = names(freqsNLsgmplus),type = "FS")
  expect_equal(nrow(x),n)
  expect_equal(ncol(x)/2,length(freqsNLsgmplus))
  expect_false(any(is.na(x)))
  
  # with theta>0
  x <- sample.relatives(x = x0,N = n,markers = names(freqsNLsgmplus),type = "FS",theta = 0.1)
  expect_equal(nrow(x),n)
  expect_equal(ncol(x)/2,length(freqsNLsgmplus))
  expect_false(any(is.na(x)))
  
  # with theta=1
  x <- sample.relatives(x = x0,N = n,markers = names(freqsNLsgmplus),type = "FS",theta = 1)
  expect_equal(nrow(x),n)
  expect_equal(ncol(x)/2,length(freqsNLsgmplus))
  expect_false(any(is.na(x)))
  
  # sampling full siblings of many x0
  x0 <- sample.profiles(N = n,freqs = freqsNLngm)
  x <- sample.relatives(x = x0,N = 1,type = "FS",markers = names(freqsNLsgmplus))
  expect_equal(nrow(x),n)
  expect_equal(ncol(x)/2,length(freqsNLsgmplus))
  expect_false(any(is.na(x)))
  
  x0 <- sample.profiles(N = n,freqs = freqsNLngm)
  x <- sample.relatives(x = x0,N = 1,type = "FS",markers = names(freqsNLsgmplus),theta = 0.1)
  expect_equal(nrow(x),n)
  expect_equal(ncol(x)/2,length(freqsNLsgmplus))
  expect_false(any(is.na(x)))  
})




