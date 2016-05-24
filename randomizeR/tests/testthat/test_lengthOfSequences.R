############################################################
# ---------------------------------------------------------#
# Test for the random Sequence Generation                  #
# Does the randomization sequence have the correct length? #                
# ---------------------------------------------------------#
############################################################
context("Correct length of randomization sequences")

test_that("randomization designs produce sequences of correct length",{
  # We test that the output matrix has the correct number of columns and thus that the
  # randomization sequence has the length set earlier (coincides with number of patients).

  # If simulation using very large numbers is desired:
  # nMax       <- c(60, 1200, 12000) 
  # seedMax    <- c(1000, 10000, 100000)
  # blockNrMax <- c(10, 100, 1000) 
  # blockMax   <- c(20, 500, 2000)
  
  nMax        <- 60
  seedMax     <- .Machine$integer.max
  blockNrMax  <- 8
  blockMax    <- 24
  
  for(i in 1:length(nMax)){
    N      <- sample(seq(2, nMax[i], 2), 1)      # Sample number of patients
    mti    <- sample(N/2, 1)                     # Sample maximum tolerated imbalance
    p      <- sample(seq(0.5001, 1, 0.05), 1)    # Sample biased coin parameter
    seed   <- sample(seedMax[i], 1)              # Sample seed
    nr     <- sample(blockNrMax[i],1)            # Sample number of blocks 
    blocks <- sample(seq(2, blockMax[i], 2), nr) # Sample blocks
    gamma  <- sample(50, 1)                      # Sample parameter for bbcd
    a      <- sample(50, 1)                      # Sample parameter for abcd
    rho    <- sample(50, 1)                      # Sample parameter for gbcd
    
    # # # # # # # # # # # # # # # # # # # 
    # Tests for K = 2
    
    # 1. Test for complete randomization
    output1 <- genSeq(crPar(N = N), seed = seed) 
    expect_equal(ncol(getRandList(output1)), N)
    
    # 2. Test for Random Allocation Rule 
    output1 <- genSeq(rarPar(N = N), seed = seed) 
    expect_equal(ncol(getRandList(output1)), N)
    
    # 3. Test for Permuted Block Randomization
    # PBR
    output1 <- genSeq(pbrPar(bc = blocks), seed = seed)
    expect_equal(ncol(getRandList(output1)), sum(blocks))
    
    #RPBR
    output2 <- genSeq(rpbrPar(rb = blocks, N = N), seed = seed)
    expect_equal(ncol(getRandList(output2)), N)

    # 4. Test for Efron's Biased Coin Desgin
    output1 <- genSeq(ebcPar(N = N, p = p), seed = seed) 
    expect_equal(ncol(getRandList(output1)), N)
    
    # 5. Test for Big Stick Design
    output1 <- genSeq(bsdPar(N = N, mti = mti), seed = seed) 
    expect_equal(ncol(getRandList(output1)), N)
    
    # 6. Test for Maximal Procedure
    output1 <- genSeq(mpPar(N = N, mti = mti), seed = seed) 
    expect_equal(ncol(getRandList(output1)), N)
    
    # 7. Test for Truncated Binomial Design
    output1 <- genSeq(tbdPar(bc = blocks), seed = seed)
    expect_equal(ncol(getRandList(output1)), sum(blocks))
    
    output2 <- genSeq(rtbdPar(N = N, rb = blocks), seed = seed)
    expect_equal(ncol(getRandList(output2)), N)
    
    # 8. Test for Urn Design
    ini <- sample(seq(2, 20, 2), 1)   # Sample initial urn composition
    add <- sample(seq(2, 20, 2), 1)   # Sample number of balls that are added to urn each step
    output1 <- genSeq(udPar(N = N, ini = ini, add = add), seed = seed)
    expect_equal(ncol(getRandList(output1)), N)
    
    # 9. Test for Hadamard Randomization
    output1 <- genSeq(hadaPar(N = N), seed = seed)
    expect_equal(ncol(getRandList(output1)), N)
    
    # 10. Test for Generalized Biased Coin Design
    output1 <- genSeq(gbcdPar(N, rho), seed = seed)
    expect_equal(ncol(getRandList(output1)), N)
    
    # 11. Test for Adjustable Biased Coin Design
    output1 <- genSeq(abcdPar(N, a), seed = seed)
    expect_equal(ncol(getRandList(output1)), N)
    
    # 12. Test for Bayesian Biased Coin Design
    output1 <- genSeq(bbcdPar(N, gamma), seed = seed)
    expect_equal(ncol(getRandList(output1)), N)
    
    # # # # # # # # # # # # # # # # # # # 
    # Test for K = 3
    
    N      <- sample(seq(3, nMax[i], 3), 1)      # Sample number of patients        
    blocks <- sample(seq(3, blockMax[i], 3), nr) # Sample blocks
    
    # ratio = c(1, 1)
    
    # 1. Test for complete randomization
    output1 <- genSeq(crPar(N = N, K = 3), seed = seed) 
    expect_equal(ncol(getRandList(output1)), N)
    
    # 2. Test for Random Allocation Rule 
    output1 <- genSeq(rarPar(N = N, K = 3), seed = seed) 
    expect_equal(ncol(getRandList(output1)), N)
    
    # 3. Test for Permuted Block Randomization
    # PBR
    output1 <- genSeq(pbrPar(bc = blocks, K = 3), seed = seed)
    expect_equal(ncol(getRandList(output1)), sum(blocks))
    
    #RPBR
    output2 <- genSeq(rpbrPar(rb = blocks, N = N, K = 3), seed = seed)
    expect_equal(ncol(getRandList(output2)), N)
    
#     # ratio != c(1, 1, 1)
#     ratio <- c(1, 2, 3)
#     N <- 24
#     
#     # 1. Test for complete randomization
#     output1 <- genSeq(crPar(N = N, K = 3, ratio = ratio), seed = seed) 
#     expect_equal(ncol(getRandList(output1)), N)
#     
#     # 2. Test for Random Allocation Rule 
#     output1 <- genSeq(rarPar(N = N, K = 3, ratio = ratio), seed = seed) 
#     expect_equal(ncol(getRandList(output1)), N)
#     
#     # 3. Test for Permuted Block Randomization
#     # PBR
#     output1 <- genSeq(pbrPar(bc = blocks, K = 3, ratio = ratio), seed = seed)
#     expect_equal(ncol(getRandList(output1)), sum(blocks))
#     
#     #RPBR
#     output2 <- genSeq(rpbrPar(rb = blocks, N = N, K = 3, ratio = ratio), seed = seed)
#     expect_equal(ncol(getRandList(output2)), N)
    
  }
})