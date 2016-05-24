
# Totally unimodular matrix, reduces to nothing

test_that("total unimodularity detection",{
  
  A <- matrix(c(
   1,1,0,0,0,
   -1,0,0,1,0,
   0,0,01,1,0,
   0,0,0,-1,1),nrow=5)
  expect_true(is_totally_unimodular(A))
  
  # Totally unimodular matrix, by Heller-Tompson criterium
  A <- matrix(c(
   1,1,0,0,
   0,0,1,1,
   1,0,1,0,
   0,1,0,1),nrow=4)
  expect_true(is_totally_unimodular(A))
  
  # Totally unimodular matrix, by Raghavachari recursive criterium
  A <- matrix(c(
      1,1,1,
      1,1,0,
      1,0,1))
  expect_true(is_totally_unimodular(A))

})


