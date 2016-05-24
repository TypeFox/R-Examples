# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# --------------------------------------------------------------- #
# Test whether frequency of treatments is approximately equal to  #
# the allocation ratio                                            #
# --------------------------------------------------------------- #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #  

context("Allocation ratio")

test_that("the frequency of treatments coincide with allocation ratio",{
  possibleN <- c(4, 6, 8)
  possibleBlocks <- list(c(4, 4), c(6, 6), c(4, 8))
  possibleRatios <- list(c(1, 1), c(1, 2), c(2,1), c(1,3))
  
  N       <- possibleN[3]                       # Number of patients
  n       <- 50000                              # Number of generated sequences
  alpha   <- 0.01                               # = 1 - confidence level
  p       <- sample(seq(0.5001, 1, 0.05), 1)    # biasesd coin parameter
  blocks  <- unlist(possibleBlocks[2])          # blocks for block based procedures
  mti     <- sample(N/2, 1)                     # Maximal Tolerated Imbalance
  ratio1  <- unlist(possibleRatios[3])          # allocation ratio for block based procedures
  r1      <- ratio1[2]/sum(ratio1)              # 'Theoretical' ratio for block based procedures
  ratio2  <- unlist(possibleRatios[4])          # allocation ratio for non-block based procedures
  r2      <- ratio2[2]/sum(ratio2)              # 'Theoretical' ratio for nonblock based procedures

  # # # # # # # # # # # # # # # # # # # #
  # 1. Test for complete randomization
  output <- genSeq(crPar(N = N, ratio = ratio2), n)
  M <- output$M
  # empirical ratio
  ratioFreq <- apply(M, 2, function(x) sum(x)/n)
  # tolerance is computed by formula from approximate standard confidence intervall
  tol <- qnorm(1-alpha/2)*sqrt(r2*(1-r2)/n)
  # difference of emprical and 'theoretical' ratio
  diff <- abs(rep(r2, N)-ratioFreq)
  # Apply the test function expect_less_than to every member of the list
  # (test whether empirical ratio is included in confidence interval)
  # lapply(diff, function(x) expect_less_than(x, tol))
  for (i in 1:length(diff)){
    expect_less_than(diff[i], tol)
  }
  
  # # # # # # # # # # # # # # # # # # # #
  # 2. Test for Random Allocation Rule 
  output <- genSeq(rarPar(N = N, ratio = ratio2), n)
  M <- output$M
  # empirical ratio
  ratioFreq <- apply(M, 2, function(x) sum(x)/n)
  # tolerance is computed by formula from approximate standard confidence intervall
  tol <- qnorm(1-alpha/2)*sqrt(r2*(1-r2)/n)
  # difference of emprical and 'theoretical' ratio
  diff <- abs(rep(r2, N)-ratioFreq)
  # Apply the test function expect_less_than to every member of the list
  # (test whether empirical ratio is included in confidence interval)
  # lapply(diff, function(x) expect_less_than(x, tol))
  for (i in 1:length(diff)){
    expect_less_than(diff[i], tol)
  }
  
  # # # # # # # # # # # # # # # # # # # # # # #
  # 3. Test for Permuted Block Randomization
  
  # PBR
  output <- genSeq(pbrPar(bc = blocks, ratio = ratio1), n)
  M <- output$M
  # empirical ratio
  ratioFreq <- apply(M, 2, function(x) sum(x)/n)
  # tolerance is computed by formula from approximate standard confidence intervall
  tol <- qnorm(1-alpha/2)*sqrt(r1*(1-r1)/n)
  # difference of emprical and 'theoretical' ratio
  diff <- abs(rep(r1, sum(blocks))-ratioFreq)
  # Apply the test function expect_less_than to every member of the list
  # (test whether empirical ratio is included in confidence interval)
  # lapply(diff, function(x) expect_less_than(x, tol))
  for (i in 1:length(diff)){
    expect_less_than(diff[i], tol)
  }
  
  # RPBR
  output <- genSeq(rpbrPar(rb = blocks, N = N, ratio = ratio1), n)
  M <- output$M
  # empirical ratio
  ratioFreq <- apply(M, 2, function(x) sum(x)/n)
  # tolerance is computed by formula from approximate standard confidence intervall
  tol <- qnorm(1-alpha/2)*sqrt(r1*(1-r1)/n)
  # difference of emprical and 'theoretical' ratio
  diff <- abs(rep(r1, N)-ratioFreq)
  # Apply the test function expect_less_than to every member of the list
  # (test whether empirical ratio is included in confidence interval)
  # lapply(diff, function(x) expect_less_than(x, tol))
  for (i in 1:length(diff)){
    expect_less_than(diff[i], tol)
  }
  
  # # # # # # # # # # # # # # # # # # # # # # 
  # 4. Test for Efron's Biased Coin Desgin
  # here always ratio = c(1,1) since EBC aims balance at the end

  output <- genSeq(ebcPar(N = N, p = p), n)
  M <- output$M
  # empirical ratio
  ratioFreq <- apply(M, 2, function(x) sum(x)/n)
  # tolerance is computed by formula from approximate standard confidence intervall
  tol <- qnorm(1-alpha/2)*sqrt(0.5*(1-0.5)/n)
  # difference of emprical and 'theoretical' ratio
  diff <- abs(rep(0.5, N)-ratioFreq)
  # Apply the test function expect_less_than to every member of the list
  # (test whether empirical ratio is included in confidence interval)
  # lapply(diff, function(x) expect_less_than(x, tol))
  for (i in 1:length(diff)){
    expect_less_than(diff[i], tol)
  }

  
  # # # # # # # # # # # # # # # # # # # # # #
  # 5. Test for Big Stick Design -> NO TEST
  
  
  # # # # # # # # # # # # # # # # #
  # 6. Test for Maximal Procedure 
  output <- genSeq(mpPar(N = N, mti = mti, ratio = ratio2), n)
  M <- output$M
  # empirical ratio
  ratioFreq <- apply(M, 2, function(x) sum(x)/n)
  # tolerance is computed by formula from approximate standard confidence intervall
  tol <- qnorm(1-alpha/2)*sqrt(r2*(1-r2)/n)
  # difference of emprical and 'theoretical' ratio
  # back to diff <- abs(rep(r, N)-ratioFreq) 
  diff <- abs(rep(r2, N)-ratioFreq)
  # Apply the test function expect_less_than to every member of the list
  # (test whether empirical ratio is included in confidence interval)
  lapply(diff, function(x) expect_less_than(x, tol))
  
  
  # # # # # # # # # # # # # # # # # # # # # 
  # 7. Test for Truncated Binomial Design
  # here only ratio = c(1, 1) since other ratios are not supported

  # TBD
  output <- genSeq(tbdPar(bc = blocks), n)
  M <- output$M
  # empirical ratio
  ratioFreq <- apply(M, 2, function(x) sum(x)/n)
  # tolerance is computed by formula from approximate standard confidence intervall
  tol <- qnorm(1-alpha/2)*sqrt(0.5*(1-0.5)/n)
  # difference of emprical and 'theoretical' ratio
  diff <- abs(rep(0.5, sum(blocks))-ratioFreq)
  # Apply the test function expect_less_than to every member of the list
  # (test whether empirical ratio is included in confidence interval)
  # lapply(diff, function(x) expect_less_than(x, tol))
  for (i in 1:length(diff)){
    expect_less_than(diff[i], tol)
  }
    
  # RTBD   
  output <- genSeq(rtbdPar(N = N, rb = blocks), n)
  M <- output$M
  # empirical ratio
  ratioFreq <- apply(M, 2, function(x) sum(x)/n)
  # tolerance is computed by formula from approximate standard confidence intervall
  tol <- qnorm(1-alpha/2)*sqrt(0.5*(1-0.5)/n)
  # difference of emprical and 'theoretical' ratio
  diff <- abs(rep(0.5, N)-ratioFreq)
  # Apply the test function expect_less_than to every member of the list
  # (test whether empirical ratio is included in confidence interval)
  # lapply(diff, function(x) expect_less_than(x, tol))
  for (i in 1:length(diff)){
    expect_less_than(diff[i], tol)
  }

  
  # # # # # # # # # # # # # # # # # # # # # 
  # 8. Test for Urn Design -> NO TEST
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # 9. Test for Hadamard Randomization => not necessary!?
  
})
