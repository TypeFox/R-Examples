library(DNAprofiles)
context("Profiles object construction")

test_that(desc = "Constructor works",{
  data(freqsNLngm,freqsNLsgmplus)
  
  # construct a profile from characters
  x1.list <- list(D1S1656 = "17.3/11", D2S441 = "11/14", D2S1338 = "19/18", 
             D3S1358 = "14/17", FGA = "26/21", D8S1179 = "12/12", D10S1248 = "15/13", 
             TH01 = "6/9.3", VWA = "18/18", D12S391 = "18/17", D16S539 = "13/13", 
             D18S51 = "12/15", D19S433 = "15/13.2", D21S11 = "29/30", 
             D22S1045 = "16/15")
  x1 <- profiles(x = x1.list,labels = get.labels(freqsNLngm))
  # and convert back to characters
  x1.chars <- profiles.to.chars(x1,freqs = freqsNLngm)
  expect_true(all(c(x1.list,recursive = TRUE)==x1.chars))

  # same for a sampled profile
  set.seed(123)
  y1 <- sample.profiles(N = 1,freqsNLngm)
  y1.chars <- profiles.to.chars(y1)
  y1.list <- sapply(colnames(y1.chars),function(nm) unname(y1.chars[,nm]),simplify = FALSE)
  y2 <- profiles(y1.list,labels = get.labels(freqsNLngm))
  
  expect_true(identical(as.vector(y1),  as.vector(y2)))
  expect_true(identical(colnames(y1), colnames(y2)  ))
  
  # same for many sampled profiles
  y1 <- sample.profiles(N = 1e2,freqsNLngm)
  y1.chars <- profiles.to.chars(y1)
  y1.list <- sapply(colnames(y1.chars),function(nm) unname(y1.chars[,nm]),simplify = FALSE)
  y2 <- profiles(y1.list,labels = get.labels(freqsNLngm))
  
  expect_true(identical(as.vector(y1),  as.vector(y2)))
  expect_true(identical(colnames(y1), colnames(y2)  ))
  
})