###############################################
# --------------------------------------------#
# Test for the random Sequence Generation     #
# Is the correct seed used?                   #
# K = 3                                       #
# --------------------------------------------#
###############################################
context("Correct Seed for K = 3")

test_that("correct random seed is used for K = 3", {
  N      <- sample(seq(3, 52, 3), 1)           # Sample number of patients
  r      <- sample(30, 1)                      # Sample number of randomisation sequences
  mti    <- sample(ceiling(N/2), 1)            # Sample maximum tolerated imbalance
  p      <- sample(seq(0.5001, 1, 0.05), 1)    # Sample biased coin parameter
  nr     <- sample(8, 1)                       # Sample number of blocks 
  blocks <- sample(seq(3, 24, 3), nr)          # Sample blocks
  seed   <- sample(.Machine$integer.max, 1)    # Sample seed
  
  # 1. Test for complete randomization
  output1 <- genSeq(crPar(N = N, K = 3), r, seed) 
  expect_equal(output1$seed, seed)
  
  # 2. Test for Random Allocation Rule 
  output1 <- genSeq(rarPar(N = N, K = 3), r, seed) 
  expect_equal(output1$seed, seed)
  
  # 3. Test for Permuted Block Randomization
  # note that N is not needed here
  output1 <- genSeq(pbrPar(bc = blocks, K = 3), r, seed)
  expect_equal(output1$seed, seed)
  
  output2 <- genSeq(rpbrPar(rb = blocks, N = N, K = 3), r, seed)
  expect_equal(output2$seed, seed)
})
