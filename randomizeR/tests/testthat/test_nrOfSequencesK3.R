###############################################################
# ------------------------------------------------------------#
# Tests for the random Sequence Generation                    #
# Does the Output have correct number of generated sequences? #
# K = 3                                                       #
# ------------------------------------------------------------#
###############################################################


context("Correct number of generated sequences for K = 3")

test_that("Output has correct number of sequences (without seed) for K = 3.", {
  # We test that the matrix of randomization sequences has the correct number of rows and 
  # thus coincides with the number of randomization sequences set earlier.
  
  N      <- sample(seq(3, 52, 3), 1)           # Sample number of patients
  r      <- sample(30, 1)                      # Sample number of randomisation sequences
  mti    <- sample(ceiling(N/2), 1)            # Sample maximum tolerated imbalance
  p      <- sample(seq(0.5001, 1, 0.05), 1)    # Sample biased coin parameter
  nr     <- sample(8, 1)                       # Sample number of blocks 
  blocks <- sample(seq(3, 24, 3), nr)          # Sample blocks
  
  # 1. Test for complete randomization
  output1 <- genSeq(crPar(N = N, K = 3), r = r) # most probably more than one sequence
  expect_equal(nrow(getRandList(output1)), r)
  
  # 2. Test for Random Allocation Rule 
  output1 <- genSeq(rarPar(N = N, K = 3), r = r) 
  expect_equal(nrow(getRandList(output1)), r)
  
  # 3. Test for Permuted Block Randomization
  # note that N is not needed here
  output1 <- genSeq(pbrPar(bc = blocks, K = 3), r = r)  # no seed
  expect_equal(nrow(getRandList(output1)), r)
  
  output2 <- genSeq(rpbrPar(rb = blocks, N = N, K = 3), r = r) 
  expect_equal(nrow(getRandList(output2)), r)
})

test_that("Output has correct number of sequences (with seed) for K = 3.", {
  # We test that the matrix of randomization sequences has the correct number of rows 
  # and thus coincides with the number randomization sequences set earlier.
  
  N      <- sample(seq(3, 52, 3), 1)           # Sample number of patients
  r      <- sample(30, 1)                      # Sample number of randomisation sequences
  mti    <- sample(ceiling(N/2), 1)            # Sample maximum tolerated imbalance
  p      <- sample(seq(0.5001, 1, 0.05), 1)    # Sample biased coin parameter
  nr     <- sample(8, 1)                       # Sample number of blocks 
  blocks <- sample(seq(3, 24, 3), nr)          # Sample blocks
  seed   <- sample(.Machine$integer.max, 1)    # Sample seed
  
  # 1. Test for complete randomization
  output1 <- genSeq(crPar(N = N, K = 3), r = r, seed = seed) 
  expect_equal(nrow(getRandList(output1)), r)
  
  # 2. Test for Random Allocation Rule 
  output1 <- genSeq(rarPar(N = N, K = 3), r = r, seed = seed) 
  expect_equal(nrow(getRandList(output1)), r)
  
  # 3. Test for Permuted Block Randomization
  # note that N is not needed here
  output1 <- genSeq(pbrPar(bc = blocks, K = 3), r = r, seed = seed)
  expect_equal(nrow(getRandList(output1)), r)
  
  output2 <- genSeq(rpbrPar(rb = blocks, N = N, K = 3), r = r, seed = seed)
  expect_equal(nrow(getRandList(output2)), r)
})
