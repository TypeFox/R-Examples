library(DNAprofiles)
context("Heterozygous/homozygous")

test_that(desc = "Heterozygous/homozygous checks",{
  # basic test
  x <- structure(t(c(1L,1L,1L,2L)),dimnames=list(NULL,c("locus1.1","locus1.2","locus2.1","locus2.2")))
  expect_equal(unname(homozygous(x)),t(c(TRUE,FALSE)))
  expect_equal(unname(heterozygous(x)),!t(c(TRUE,FALSE)))
  
  data(freqsNLsgmplus)
  set.seed(100)
  for(n in c(0,1,10,100)){
    x <- sample.profiles(N = n,freqs = freqsNLsgmplus)
    x.het <- heterozygous(x)
    x.hom <- homozygous(x)
    
    expect_equal(TRUE,all(x.het|x.hom))
    expect_equal(dim(x),dim(x.het)*c(1,2))  
    expect_equal(dim(x),dim(x.hom)*c(1,2))
    expect_equal(TRUE,all((rowSums(x.het)+rowSums(x.hom))==length(freqsNLsgmplus)))      
  }
})