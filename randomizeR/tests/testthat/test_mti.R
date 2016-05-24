# # # # # # # # # # # # # # # # # # # # # # # # # # #
# ------------------------------------------------- #
# Test for Maximal Procedure and Big Stick Design   #
# Do the methods really comply with MTI?            #
# ------------------------------------------------- #
# # # # # # # # # # # # # # # # # # # # # # # # # # # 


context("Comply with MTI")

test_that("the randomization methods comply with the prescribed MTI",{
  N <- sample(seq(2, 20, 2), 1)
  mti <- sample(N/2, 1)
  p <- sample(seq(0.5, 1, 0.01), 1)
  nr <- sample(10, 1)
  
  # Maximal Procedure
  mpSeq <- genSeq(mpPar(N = N, mti = mti), nr)
  # Produce random walk
  randomWalk <- t(apply(mpSeq$M, 1, function(x) { cumsum(c(2*x-1)) }))
  # Now check, whether each entry does not exceed the mti
  expect_true(all(abs(randomWalk) <= mti))
  
  
  # Big Stick Design
  bsdSeq <- genSeq(bsdPar(N = N, mti = mti), nr)
  randomWalk <- t(apply(bsdSeq$M, 1, function(x) { cumsum(c(2*x-1)) }))
  # Now check, whether each entry does not exceed the mti
  expect_true(all(abs(randomWalk) <=  mti))
  
  # Chen's Design
  chenSeq <- genSeq(chenPar(N = N, mti = mti, p = p))
  randomWalk <- t(apply(chenSeq$M, 1, function(x) { cumsum(c(2*x-1)) }))
  expect_true(all(abs(randomWalk) <= mti))
})