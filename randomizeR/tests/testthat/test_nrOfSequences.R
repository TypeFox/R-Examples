###############################################################
# ------------------------------------------------------------#
# Tests for the random Sequence Generation                    #
# Does the Output have correct number of generated sequences? #
# ------------------------------------------------------------#
###############################################################


context("Correct number of generated sequences")

test_that("Output has correct number of sequences (without seed).", {
  # We test that the matrix of randomization sequences has the correct number of rows and 
  # thus coincides with the number of randomization sequences set earlier.
  
  N      <- sample(seq(2, 50, 2), 1)           # Sample number of patients
  r      <- sample(30, 1)                      # Sample number of randomisation sequences
  mti    <- sample(N/2, 1)                     # Sample maximum tolerated imbalance
  p      <- sample(seq(0.5001, 1, 0.05), 1)    # Sample biased coin parameter
  nr     <- sample(10,1)                       # Sample number of blocks 
  blocks <- sample(seq(2, 20, 2), nr)          # Sample blocks
  gamma  <- sample(50, 1)                      # Sample parameter for bbcd
  a      <- sample(50, 1)                      # Sample parameter for abcd
  rho    <- sample(50, 1)                      # Sample parameter for gbcd
  
  # 1. Test for complete randomization
  output1 <- genSeq(crPar(N = N), r = r) # most probably more than one sequence
  expect_equal(nrow(getRandList(output1)), r)
  
  # 2. Test for Random Allocation Rule 
  output1 <- genSeq(rarPar(N = N), r = r) 
  expect_equal(nrow(getRandList(output1)), r)
  
  # 3. Test for Permuted Block Randomization
  # note that N is not needed here
  output1 <- genSeq(pbrPar(bc = blocks), r = r)  # no seed
  expect_equal(nrow(getRandList(output1)), r)
  
  output2 <- genSeq(rpbrPar(rb = blocks, N = N), r = r) 
  expect_equal(nrow(getRandList(output2)), r)
  
  # 4. Test for Efron's Biased Coin Desgin
  output1 <- genSeq(ebcPar(N = N, p = p), r = r) 
  expect_equal(nrow(getRandList(output1)), r)
  
  # 5. Test for Big Stick Design
  output1 <- genSeq(bsdPar(N = N, mti = mti), r = r) 
  expect_equal(nrow(getRandList(output1)), r)
  
  # 6. Test for Maximal Procedure
  output1 <- genSeq(mpPar(N = N, mti = mti), r = r) 
  expect_equal(nrow(getRandList(output1)), r)
  
  # 7. Test for Truncated Binomial Design
  output1 <- genSeq(tbdPar(bc = blocks), r = r) 
  expect_equal(nrow(getRandList(output1)), r)
  
  output2 <- genSeq(rtbdPar(N = N, rb = blocks), r = r)
  expect_equal(nrow(getRandList(output2)), r)
  
  # 8. Test for Urn Design
  ini <- sample(seq(2, 20, 2), 1)   # Sample initial urn composition
  add <- sample(seq(2, 20, 2), 1)   # Sample number of balls that are added to urn each step
  output1 <- genSeq(udPar(N = N, ini = ini, add = add), r = r) # no seed
  expect_equal(nrow(getRandList(output1)), r)
  
  # 9. Test for Hadamard Randomization
  output1 <- genSeq(hadaPar(N = N), r = r) 
  expect_equal(nrow(getRandList(output1)), r)
  
  # 10. Test for Generalized Biased Coin Design
  output1 <- genSeq(gbcdPar(N, rho), r = r)
  expect_equal(nrow(getRandList(output1)), r)
  
  # 11. Test for Adjustable Biased Coin Design
  output1 <- genSeq(abcdPar(N, a), r = r)
  expect_equal(nrow(getRandList(output1)), r)
  
  # 12. Test for Bayesian Biased Coin Design
  output1 <- genSeq(bbcdPar(N, gamma), r = r)
  expect_equal(nrow(getRandList(output1)), r)
})

test_that("Output has correct number of sequences (with seed)", {
  # We test that the matrix of randomization sequences has the correct number of rows 
  # and thus coincides with the number randomization sequences set earlier.
  
  N      <- sample(seq(2, 50, 2), 1)           # Sample number of patients
  r      <- sample(30, 1)                      # Sample number of randomisation sequences
  mti    <- sample(N/2, 1)                     # Sample maximum tolerated imbalance
  p      <- sample(seq(0.5001, 1, 0.05), 1)    # biased coin parameter
  seed   <- sample(.Machine$integer.max, 1)    # Sample seed
  nr     <- sample(10,1)                       # sample number of blocks 
  blocks <- sample(seq(2, 20, 2), nr)          # sample blocks
  gamma  <- sample(50, 1)                      # Sample parameter for bbcd
  a      <- sample(50, 1)                      # Sample parameter for abcd
  rho    <- sample(50, 1)                      # Sample parameter for gbcd
  
  # 1. Test for complete randomization
  output1 <- genSeq(crPar(N = N), r = r, seed = seed) 
  expect_equal(nrow(getRandList(output1)), r)
  
  # 2. Test for Random Allocation Rule 
  output1 <- genSeq(rarPar(N = N), r = r, seed = seed) 
  expect_equal(nrow(getRandList(output1)), r)
  
  # 3. Test for Permuted Block Randomization
  # note that N is not needed here
  output1 <- genSeq(pbrPar(bc = blocks), r = r, seed = seed)
  expect_equal(nrow(getRandList(output1)), r)
  
  output2 <- genSeq(rpbrPar(rb = blocks, N = N), r = r, seed = seed)
  expect_equal(nrow(getRandList(output2)), r)
  
  # 4. Test for Efron's Biased Coin Desgin
  output1 <- genSeq(ebcPar(N = N, p = p), r = r, seed = seed) 
  expect_equal(nrow(getRandList(output1)), r)
  
  # 5. Test for Big Stick Design
  output1 <- genSeq(bsdPar(N = N, mti = mti), r = r, seed = seed) 
  expect_equal(nrow(getRandList(output1)), r)
  
  # 6. Test for Maximal Procedure
  output1 <- genSeq(mpPar(N = N, mti = mti), r = r, seed = seed) 
  expect_equal(nrow(getRandList(output1)), r)
  
  # 7. Test for Truncated Binomial Design
  output1 <- genSeq(tbdPar(bc = N), r = r, seed = seed)
  expect_equal(nrow(getRandList(output1)), r)
  
  output2 <- genSeq(rtbdPar(N = N, rb = blocks), r = r, seed = seed)
  expect_equal(nrow(getRandList(output2)), r)
  
  # 8. Test for Urn Design
  ini <- sample(seq(2, 20, 2), 1)   # Sample initial urn composition
  add <- sample(seq(2, 20, 2), 1)   # Sample number of balls that are added to urn each step
  output1 <- genSeq(udPar(N = N, ini = ini, add = add), r = r, seed = seed)
  expect_equal(nrow(getRandList(output1)), r)
  
  # 9. Test for Hadamard Randomization
  output1 <- genSeq(hadaPar(N = N), r = r, seed = seed)
  expect_equal(nrow(getRandList(output1)), r)
  
  # 10. Test for Generalized Biased Coin Design
  output1 <- genSeq(gbcdPar(N, rho), r = r, seed = seed)
  expect_equal(nrow(getRandList(output1)), r)
  
  # 11. Test for Adjustable Biased Coin Design
  output1 <- genSeq(abcdPar(N, a), r = r, seed = seed)
  expect_equal(nrow(getRandList(output1)), r)
  
  # 12. Test for Bayesian Biased Coin Design
  output1 <- genSeq(bbcdPar(N, gamma), r = r, seed = seed)
  expect_equal(nrow(getRandList(output1)), r)
})

