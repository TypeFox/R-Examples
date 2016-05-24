# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ----------------------------------------------------------- #
# Test for block based randomization methods                  #
# that there is the desired balance balance within each block #                                    
# ----------------------------------------------------------- #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

context("Balance within blocks")

test_that("there is an approximate balance Within blocks",{
  N      <- sample(seq(2, 50, 2), 1)     # Sample number of patients
  nr     <- sample(10, 1)                # Sample number of blocks
  blocks <- sample(seq(4, 22, 2), nr)    # Sample blocks
  
                            
  # # # # # # # # # 
  # WITHout RATIO
  # # # # # # # # #
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # Permuted Block Randomization WITHOUT ratio, K = 2
  pbrSeq   <- genSeq(pbrPar(bc = blocks))
  nrBlocks <- length(pbrSeq$bc)
  j <-  1
  for(i in 1:nrBlocks){
    iBlockLength <- pbrSeq$bc[i]
    # Here: test for balance
    expect_equal(sum(pbrSeq$M[j:(j+iBlockLength-1)]), iBlockLength/2)
    j <- j + iBlockLength
  }
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # Randomized Permuted Block Randomization WITHOUT ratio, K = 2

  rpbrSeq <- genSeq(rpbrPar(rb = blocks, N = N))
  nrBlocks <- length(rpbrSeq$bc[[1]])
  
  j <-  1
  if(N > rpbrSeq$bc[[1]][1]){
    for(i in 1:nrBlocks){
      iBlockLength <- rpbrSeq$bc[[1]][i]
      # Only check, if  full blocks are used, i.e. the blocks completely filled
      if(j + iBlockLength - 1 <= N){
        # Here: test for balance
        expect_equal(sum(rpbrSeq$M[j:(j + iBlockLength - 1)]), iBlockLength/2)
        j <- j + iBlockLength
      }
    }
  }
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # Truncated Binomial Design WITHOUT ratio, K = 2
  tbdSeq <- genSeq(tbdPar(bc = blocks))
  nrBlocks <- length(tbdSeq$bc)
  
  j <-  1
  for(i in 1:nrBlocks){
    iBlockLength <- tbdSeq$bc[i]
      # Here: test for balance
      expect_equal(sum(tbdSeq$M[j:(j + iBlockLength - 1)]), iBlockLength/2)
      j <- j + iBlockLength
  }
  
  rtbdSeq <- genSeq(rtbdPar(N = N, rb = blocks))
  nrBlocks <- length(rtbdSeq$bc[[1]])
  
  j <-  1
  if(N > rtbdSeq$bc[[1]][1]){
    for(i in 1:nrBlocks){
      iBlockLength <- rtbdSeq$bc[[1]][i]
      # Only check, if  full blocks are used, i.e. the blocks completely filled
      if(j + iBlockLength - 1 <= N){
        # Here: test for balance
        expect_equal(sum(rtbdSeq$M[j:(j + iBlockLength - 1)]), iBlockLength/2)
        j <- j + iBlockLength
      }
    }
  }
  
  # # # # # # # # # 
  # WITH RATIO
  # # # # # # # # #
  
  # # # # # # # # # # # # # # # # # # # # # # # # # 
  # Permuted Block Randomization WITH ratio, K = 2
  
  # einfaches Sampling von ratio:
  # Test for ratio = c(1,3), c(2,4), c(3,5)
  possibleRatios <- list(c(1,3), c(2,4), c(5,3))
  ratio          <- unlist(sample(possibleRatios, 1))
  
  # Sampling der Bloecke 
  nrBlocks       <- sample(5, 1)
  allBlocks      <- unlist(lapply(1:8, function(x) {sum(ratio)*x}))
  blocks         <- sample(allBlocks, nrBlocks)
  
  # Generates randomization sequence pbr
  pbrSeq <- genSeq(pbrPar(bc = blocks, ratio = ratio))

  j <-  1
  for(i in 1:nrBlocks){
    iBlockLength <- pbrSeq$bc[i]
    factor <- iBlockLength/sum(ratio)
    # Here: test for balance (with respect to ratio)
    expect_equal(sum(pbrSeq$M[j:(j+iBlockLength-1)]), ratio[2]*factor)
    j <- j + iBlockLength
  }
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # Randomized Permuted Block Randomization WITH ratio, K = 2
  
  # Only check, if  full blocks are used, i.e. the blocks completely filled
  N <- sample(unlist(lapply(3:10, function(x) {sum(ratio)*x})), 1)
  rpbrSeq <- genSeq(rpbrPar(rb = blocks, N = N, ratio = ratio))
  nrBlocks <- length(rpbrSeq$bc[[1]])

  j <- 1
  if(N > rpbrSeq$bc[[1]][1]){
    for(i in 1:nrBlocks){
      iBlockLength <- rpbrSeq$bc[[1]][i]
      if(j + iBlockLength - 1 <= N){
        factor <- iBlockLength/sum(ratio)
        # Here: test for balance (with respect to ratio)
        expect_equal(sum(rpbrSeq$M[j:(j+iBlockLength-1)]), ratio[2]*factor)
        j <- j + iBlockLength
      }
    }
  }
  

  # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
  # Truncated Binomial Design WITH ratio, K = 2 NOT SUPPORTED
  
})





  