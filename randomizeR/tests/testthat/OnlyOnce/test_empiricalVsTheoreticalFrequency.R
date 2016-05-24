# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# --------------------------------------------------------------- #
# Test whether empirical and and theoretical frequency coincide   #                                    #
# --------------------------------------------------------------  #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

context("Emprical vs theoretical frequency")

test_that("the empirical and theoretical frequency coincide",{
  N       <- 4                                      # Number of patients
  n       <- 100000                                 # Number of generated sequences
  alpha   <- 0.01                                   # = 1 - confidence level
  p       <- sample(seq(0.5001, 1, 0.05), 1)        # biasesd coin parameter
  blocks  <- c(4,4)                                 # blocks for block based procedures
  mti     <- sample(N/2, 1)                         # Maximal Tolerated Imbalance
  gamma   <- sample(50, 1)                          # Sample parameter for bbcd
  a       <- sample(50, 1)                          # Sample parameter for abcd
  rho     <- sample(50, 1)                          # Sample parameter for gbcd

  
  # # # # # # # # # # # # # # # # # # # #
  # 1. Test for complete randomization
  
  # Random sequence is generated
  output1 <- genSeq(crPar(N = N))
  M1 <- output1$M[1,] 
  # the corresponding theoretical frequency is computes
  p1 <- getProb(output1)
  
  # n sequences are generated
  outputN <- genSeq(crPar(N = N), n)
  Mn <- outputN$M
  # empirical frequency of sequence of above is computed
  # check which rows of matrix Mn coincide with random sequence M1, count them and 
  # compute the relative frequency
  freq <- length(which(apply(Mn, 1, function(x) all(x == M1))))/n
  # tolerance is computed by formula from approximate standard confidence intervall
  tol <- qnorm(1-alpha/2)*sqrt(p1*(1-p1)/n)   
  # check if empirical frequency is included confidence interval                                    
  expect_less_than(abs(p1-freq), tol)
  
  
  # # # # # # # # # # # # # # # # # # # #
  # 2. Test for Random Allocation Rule 
  
  # Random sequence is generated
  output1 <- genSeq(rarPar(N = N)) 
  M1 <- output1$M[1,] 
  # the corresponding theoretical frequency is computes
  p1 <- getProb(output1)
  
  # n sequences are generated
  outputN <- genSeq(rarPar(N = N), n)
  Mn <- outputN$M
  # empirical frequency of sequence of above is computed
  freq <- length(which(apply(Mn, 1, function(x) all(x == M1))))/n
  # tolerance is computed by formula from approximate standard confidence intervall
  tol <- qnorm(1-alpha/2)*sqrt(p1*(1-p1)/n)   
  # check if empirical frequency is included confidence interval                                    
  expect_less_than(abs(p1-freq), tol)
  
  
  # # # # # # # # # # # # # # # # # # # # # # #
  # 3. Test for Permuted Block Randomization
  
  # Random sequence is generated
  output1 <- genSeq(pbrPar(bc = blocks))  # no seed
  M1 <- output1$M[1,] 
  # the corresponding theoretical frequency is computes
  p1 <- getProb(output1)
  
  # n sequences are generated
  outputN <- genSeq(pbrPar(bc = blocks), n)
  Mn <- outputN$M
  # empirical frequency of sequence of above is computed
  freq <- length(which(apply(Mn, 1, function(x) all(x == M1))))/n
  # tolerance is computed by formula from approximate standard confidence intervall
  tol <- qnorm(1-alpha/2)*sqrt(p1*(1-p1)/n)   
  # check if empirical frequency is included confidence interval                                    
  expect_less_than(abs(p1-freq), tol)
  
  
  # # # # # # # # # # # # # # # # # # # # # # 
  # 4. Test for Efron's Biased Coin Desgin
  
  # Random sequence is generated
  output1 <- genSeq(ebcPar(N = N, p = p)) 
  M1 <- output1$M[1,] 
  # the corresponding theoretical frequency is computes
  p1 <- getProb(output1)
  
  # n sequences are generated
  outputN <- genSeq(ebcPar(N = N, p = p), n)
  Mn <- outputN$M
  # empirical frequency of sequence of above is computed
  freq <- length(which(apply(Mn, 1, function(x) all(x == M1))))/n
  # tolerance is computed by formula from approximate standard confidence intervall
  tol <- qnorm(1-alpha/2)*sqrt(p1*(1-p1)/n)   
  # check if empirical frequency is included confidence interval                                    
  expect_less_than(abs(p1-freq), tol)
  
  
  # # # # # # # # # # # # # # # # #
  # 5. Test for Big Stick Design
  
  # Random sequence is generated
  output1 <- genSeq(bsdPar(N = N, mti = mti)) 
  M1 <- output1$M[1,] 
  # the corresponding theoretical frequency is computes
  p1 <- getProb(output1)
  
  # n sequences are generated
  outputN <- genSeq(bsdPar(N = N, mti = mti), n)
  Mn <- outputN$M
  # empirical frequency of sequence of above is computed
  freq <- length(which(apply(Mn, 1, function(x) all(x == M1))))/n
  # tolerance is computed by formula from approximate standard confidence intervall
  tol <- qnorm(1-alpha/2)*sqrt(p1*(1-p1)/n)   
  # check if empirical frequency is included confidence interval                                    
  expect_less_than(abs(p1-freq), tol)
  
  
  # # # # # # # # # # # # # # # # #
  # 6. Test for Maximal Procedure
  
  # Random sequence is generated
  output1 <- genSeq(mpPar(N = N, mti = mti)) 
  M1 <- output1$M[1,] 
  # the corresponding theoretical frequency is computes
  p1 <- getProb(output1)
  
  # n sequences are generated
  outputN <- genSeq(mpPar(N = N, mti = mti), n)
  Mn <- outputN$M
  # empirical frequency of sequence of above is computed
  freq <- length(which(apply(Mn, 1, function(x) all(x == M1))))/n
  # tolerance is computed by formula from approximate standard confidence intervall
  tol <- qnorm(1-alpha/2)*sqrt(p1*(1-p1)/n)   
  # check if empirical frequency is included confidence interval                                    
  expect_less_than(abs(p1-freq), tol)
  
  
  # # # # # # # # # # # # # # # # # # # # # 
  # 7. Test for Truncated Binomial Design
  
  # Random sequence is generated
  output1 <- genSeq(tbdPar(bc = blocks)) 
  M1 <- output1$M[1,] 
  # the corresponding theoretical frequency is computes
  p1 <- getProb(output1)
  
  # n sequences are generated
  outputN <- genSeq(tbdPar(bc = blocks), n)
  Mn <- outputN$M
  # empirical frequency of sequence of above is computed
  freq <- length(which(apply(Mn, 1, function(x) all(x == M1))))/n
  # tolerance is computed by formula from approximate standard confidence intervall
  tol <- qnorm(1-alpha/2)*sqrt(p1*(1-p1)/n)   
  # check if empirical frequency is included confidence interval                                    
  expect_less_than(abs(p1-freq), tol)

  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # 9. Test for Hadamard Randomization => not necessary!?
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # 10. Test for Generalized Biased Coin Design
  output1 <- genSeq(gbcdPar(N, rho))
  M1 <- output1$M[1,] 
  # the corresponding theoretical frequency is computes
  p1 <- getProb(output1)
  
  # n sequences are generated
  outputN <- genSeq(gbcdPar(N, rho), n)
  Mn <- outputN$M
  # empirical frequency of sequence of above is computed
  freq <- length(which(apply(Mn, 1, function(x) all(x == M1))))/n
  # tolerance is computed by formula from approximate standard confidence intervall
  tol <- qnorm(1-alpha/2)*sqrt(p1*(1-p1)/n)   
  # check if empirical frequency is included confidence interval                                    
  expect_less_than(abs(p1-freq), tol)
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # 11. Test for Adjustable Biased Coin Design
  output1 <- genSeq(abcdPar(N, a))
  M1 <- output1$M[1,] 
  # the corresponding theoretical frequency is computes
  p1 <- getProb(output1)
  
  # n sequences are generated
  outputN <- genSeq(abcdPar(N, a), n)
  Mn <- outputN$M
  # empirical frequency of sequence of above is computed
  freq <- length(which(apply(Mn, 1, function(x) all(x == M1))))/n
  # tolerance is computed by formula from approximate standard confidence intervall
  tol <- qnorm(1-alpha/2)*sqrt(p1*(1-p1)/n)   
  # check if empirical frequency is included confidence interval                                    
  expect_less_than(abs(p1-freq), tol)
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # 12. Test for Bayesian Biased Coin Design
  output1 <- genSeq(bbcdPar(N, gamma))
  M1 <- output1$M[1,] 
  # the corresponding theoretical frequency is computes
  p1 <- getProb(output1)
  
  # n sequences are generated
  outputN <- genSeq(bbcdPar(N, gamma), n)
  Mn <- outputN$M
  # empirical frequency of sequence of above is computed
  freq <- length(which(apply(Mn, 1, function(x) all(x == M1))))/n
  # tolerance is computed by formula from approximate standard confidence intervall
  tol <- qnorm(1-alpha/2)*sqrt(p1*(1-p1)/n)   
  # check if empirical frequency is included confidence interval                                    
  expect_less_than(abs(p1-freq), tol)
})
