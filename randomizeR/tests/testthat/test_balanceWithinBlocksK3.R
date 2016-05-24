# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ----------------------------------------------------------- #
# Test for block based randomization methods                  #
# that there is the desired balance balance within each block #
# K = 3                                                       #
# ----------------------------------------------------------- #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

context("Balance within blocks for K = 3")

test_that("there is an approximate balance Within blocks for K = 3",{
  N      <- sample(seq(3, 52, 3), 1)     # Sample number of patients
  nr     <- sample(8, 1)                 # Sample number of blocks
  blocks <- sample(seq(3, 24, 3), nr)    # Sample blocks
  
  # # # # # # # # # 
  # WITHout RATIO
  # # # # # # # # #
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # Permuted Block Randomization WITHOUT ratio, K = 3
  pbrSeq   <- genSeq(pbrPar(bc = blocks, K = 3))
  nrBlocks <- length(pbrSeq$bc)
  j <-  1
  for(i in 1:nrBlocks){
    iBlockLength <- pbrSeq$bc[i]
    # Here: test for balance
    # Note: test for == 2 is not necessary, since it follows from the other tests
    expect_equal(sum(pbrSeq$M[j:(j+iBlockLength-1)]==0), iBlockLength/3)
    expect_equal(sum(pbrSeq$M[j:(j+iBlockLength-1)]==1), iBlockLength/3)
    j <- j + iBlockLength
  }
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # Randomized Permuted Block Randomization WITHOUT ratio, K = 3
  
  rpbrSeq <- genSeq(rpbrPar(rb = blocks, N = N, K = 3))
  nrBlocks <- length(rpbrSeq$bc[[1]])
  
  j <-  1
  if(N > rpbrSeq$bc[[1]][1]){
    for(i in 1:nrBlocks){
      iBlockLength <- rpbrSeq$bc[[1]][i]
      # Only check, if  full blocks are used, i.e. the blocks completely filled
      if(j + iBlockLength - 1 <= N){
        # Here: test for balance
        # Note: test for == 2 is not necessary, since it follows from the other tests
        expect_equal(sum(rpbrSeq$M[j:(j + iBlockLength - 1)] == 0), 
                     iBlockLength/3)
        expect_equal(sum(rpbrSeq$M[j:(j + iBlockLength - 1)] == 1), 
                     iBlockLength/3)
        j <- j + iBlockLength
      }
    }
  }
  
  # # # # # # # # # 
  # WITH RATIO
  # # # # # # # # #
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # Permuted Block Randomization WITH ratio, K = 3
  
  # simple sampling of ratio:
  # Test for ratio = c(1,2,3), c(5,3,1), c(7,1,4)
  possibleRatios <- list(c(1, 2, 3), c(5, 3, 1), c(7, 1, 4))
  ratio          <- unlist(sample(possibleRatios, 1))
  
  # Sampling of blocks:
  nrBlocks       <- sample(5, 1)
  allBlocks      <- unlist(lapply(1:8, function(x) {sum(ratio)*x}))
  blocks         <- sample(allBlocks, nrBlocks)
  
  # Generates randomization sequence pbr
  pbrSeq <- genSeq(pbrPar(bc = blocks, K = 3, ratio = ratio))
  
  j <-  1
  for(i in 1:nrBlocks){
    iBlockLength <- pbrSeq$bc[i]
    factor <- iBlockLength/sum(ratio)
    # Here: test for balance (with respect to ratio)
    # Note: test for == 2 is not necessary, since it follows from the other tests
    expect_equal(sum(pbrSeq$M[j:(j+iBlockLength-1)] == 0), ratio[1]*factor)
    expect_equal(sum(pbrSeq$M[j:(j+iBlockLength-1)] == 1), ratio[2]*factor)
    j <- j + iBlockLength
  }
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # Randomized Permuted Block Randomization WITH ratio, K = 3
  
  # Only check, if  full blocks are used, i.e. the blocks completely filled
  N <- sample(unlist(lapply(1:8, function(x) {sum(ratio)*x})), 1)
  rpbrSeq <- genSeq(rpbrPar(rb = blocks, N = N, K = 3, ratio = ratio))
  nrBlocks <- length(rpbrSeq$bc[[1]])
  
  j <- 1
  if (N > rpbrSeq$bc[[1]][1]) {
    for (i in 1:nrBlocks) {
      iBlockLength <- rpbrSeq$bc[[1]][i]
      if (j + iBlockLength - 1 <= N) {
        factor <- iBlockLength/sum(ratio)
        # Here: test for balance (with respect to ratio)
        # Note: test for == 2 is not necessary, since it follows from the other tests
        expect_equal(sum(rpbrSeq$M[j:(j+iBlockLength-1)] == 0), ratio[1]*factor)
        expect_equal(sum(rpbrSeq$M[j:(j+iBlockLength-1)] == 1), ratio[2]*factor)
        j <- j + iBlockLength
      }
    }
  }
})





