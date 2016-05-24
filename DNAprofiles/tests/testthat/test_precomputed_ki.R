library(DNAprofiles)
context("Precomputing KIs")

test_that(desc = "Correctness of ki.db with precomputed KIs",{
  tol <- 1e-6
  
  data(freqsNLngm)
  set.seed(123)
  
  x1 <- sample.profiles(N = 1e2,freqs = freqsNLngm)
  x2 <- sample.profiles(N = 1e3,freqs = freqsNLngm)

  M1 <- sample(names(freqsNLngm),size = 4)
  M2 <- sample(M1,size = 3)
  
  precomp <- ki.precompute(type = "FS",freqs = freqsNLngm,markers = M1)
  
  R1 <- ki.db(x = x1,db = x2,hyp.1 = "FS",markers = M2)
  R2 <- ki.db(x = x1,db = x2,markers = M2, precomputed.kis = precomp)
  
  expect_true(all.equal(R1,R2, tolerance = tol,scale = R1),TRUE)  
})
  
