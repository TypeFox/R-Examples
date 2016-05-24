###############################################
# --------------------------------------------#
# Test for the random Sequence Generation     #
# Is the correct seed used?                   #
# --------------------------------------------#
###############################################
context("Correct Seed")

test_that("correct random seed is used", {
  N      <- sample(seq(2, 50, 2), 1)           # Sample number of patients
  r      <- sample(30, 1)                      # Sample number of randomisation sequences
  mti    <- sample(N/2, 1)                     # Sample maximum tolerated imbalance
  p      <- sample(seq(0.5001, 1, 0.05), 1)    # biased coin parameter
  seed   <- sample(1000, 1)                    # Sample seed
  nr     <- sample(10,1)                       # sample number of blocks 
  blocks <- sample(seq(2, 20, 2), nr)          # sample blocks
  gamma  <- sample(50, 1)                      # Sample parameter for bbcd
  a      <- sample(50, 1)                      # Sample parameter for abcd
  rho    <- sample(50, 1)                      # Sample parameter for gbcd
  
  # 1. Test for complete randomization
  output1 <- genSeq(crPar(N), r, seed) 
  expect_equal(output1$seed, seed)
  
  # 2. Test for Random Allocation Rule 
  output1 <- genSeq(rarPar(N), r, seed) 
  expect_equal(output1$seed, seed)
  
  # 3. Test for Permuted Block Randomization
  # note that N is not needed here
  output1 <- genSeq(pbrPar(bc = blocks), r, seed)
  expect_equal(output1$seed, seed)
  
  output2 <- genSeq(rpbrPar(rb = blocks, N), r, seed)
  expect_equal(output2$seed, seed)
  
  # 4. Test for Efron's Biased Coin Desgin
  output1 <- genSeq(ebcPar(N, p), r, seed) 
  expect_equal(output1$seed, seed)
  
  # 5. Test for Big Stick Design
  output1 <- genSeq(bsdPar(N, mti), r, seed) 
  expect_equal(output1$seed, seed)
  
  # 6. Test for Maximal Procedure
  output1 <- genSeq(mpPar(N, mti), r, seed) 
  expect_equal(output1$seed, seed)
  
  # 7. Test for Truncated Binomial Design
  output1 <- genSeq(tbdPar(bc = blocks), r, seed)
  expect_equal(output1$seed, seed)

  output2 <- genSeq(rtbdPar(N, rb = blocks), r, seed)
  expect_equal(output2$seed, seed)

  # 8. Test for Urn Design
  ini <- sample(seq(2, 20, 2), 1)   # Sample initial urn composition
  add <- sample(seq(2, 20, 2), 1)   # Sample number of balls that are added to urn each step
  output1 <- genSeq(udPar(N, ini, add), r, seed)
  expect_equal(output1$seed, seed)
  
  # 9. Test for Hadamard Randomization
  output1 <- genSeq(hadaPar(N), r, seed)
  expect_equal(output1$seed, seed)
  
  # 10. Test for Generalized Biased Coin Design
  output1 <- genSeq(gbcdPar(N, rho), r, seed)
  expect_equal(output1$seed, seed)
  
  # 11. Test for Adjustable Biased Coin Design
  output1 <- genSeq(abcdPar(N, a), r, seed)
  expect_equal(output1$seed, seed)
  
  # 12. Test for Bayesian Biased Coin Design
  output1 <- genSeq(bbcdPar(N, gamma), r, seed)
  expect_equal(output1$seed, seed)
})
