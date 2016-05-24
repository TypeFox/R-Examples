# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# ----------------------------------------------------------------------------- #
# Test for Random Allocation Rule and Truncated Binomial Design                 #
# that there is the desired balance at the end of each randomization,           #
# i.e. (if no ratio given) there should be n/2 patients in each treatment group #
# (implemented for K = 2, 3)                                                    #
# ----------------------------------------------------------------------------- #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

context("Balance within treatment groups")

test_that("there is an approximate balance at the end of the randomization process",{
  
  N <- sample(seq(2, 20, 2), 1)   # Number of patients
  n1 <- sample(N-1, 1)            # sample ratio
  n2 <- N - n1                
  ratio <- c(n1, n2)          
  nr <- sample(20, 1)             # number of randomization sequences (rows of matrix)
  
  # # # # # # # # # # # # # # # # # # # # # # # #
  # Random Allocation Rule without ratio, K = 2
  rarSeq <- genSeq(rarPar(N = N), r = nr)
  # Test that the sum of each row of the matrix is equal to N/2
  expect_equal(apply(rarSeq$M, 1, function(x) sum(x)), rep(N/2, nrow(rarSeq$M)))
  
  # # # # # # # # # # # # # # # # # # # # # # # #
  # Random Allocation Rule with ratio, K = 2
  rarSeq <- genSeq(rarPar(N = N, ratio = ratio), r = nr)
  # Test that the sum of each row of the matrix is equal to the corresponding ratio
  expect_equal(apply(rarSeq$M, 1, function(x) sum(x)), rep(n2, nrow(rarSeq$M)))
  
  
  # # # # # # # # # # # # # # # # # # # # # # # # # # # #
  # Truncated Binomial Design without ratio
  tbdSeq <- genSeq(tbdPar(bc = N), r = nr)
  expect_equal(apply(tbdSeq$M, 1, function(x) sum(x)), rep(N/2, nrow(tbdSeq$M)))
  
  
  # No test for RTBD as balance can only occur if all blocks are filled and it is
  # very unlikely that all blocks are filled
  
  # # # # # #
  # K = 3
  # # # # # #
  N <- sample(seq(3, 24, 3), 1)   # Number of patients
  n1 <- sample(seq(1,N-2), 1)            # sample ratio
  n2 <- sample(seq(1,N-n1-1), 1)
  n3 <- N - n1 - n2
  ratio <- c(n1, n2, n3)          
  nr <- sample(20, 1)             # number of randomization sequences (rows of matrix)
  
  # # # # # # # # # # # # # # # # # # # # # # # #
  # Random Allocation Rule without ratio, K = 3
  rarSeq <- genSeq(rarPar(N = N, K = 3), r= nr)
  # Test that in each row of the matrix the number of occurences of 0 or 1 is equal to N/3
  # respectively, note: the check for == 2 is not necessary as it follows from the others
  expect_equal(apply(rarSeq$M, 1, function(x) sum(x == 0)), rep(N/3, nr))
  expect_equal(apply(rarSeq$M, 1, function(x) sum(x == 1)), rep(N/3, nr))
  
  # # # # # # # # # # # # # # # # # # # # # # # #
  # Random Allocation Rule with ratio, K = 3
  rarSeq <- genSeq(rarPar(N = N, K = 3, ratio = ratio), r = nr)
  # Test that in each row of the matrix the number of occurences of 0 or 1 is equal to 
  # the corresponding ratio, note: the check for == 2 is not necessary as it follows 
  # from the others
  expect_equal(apply(rarSeq$M, 1, function(x) sum(x == 0)), rep(n1, nr))
  expect_equal(apply(rarSeq$M, 1, function(x) sum(x == 1)), rep(n2, nr))
  
})