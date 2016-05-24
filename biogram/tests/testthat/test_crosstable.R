context("Crosstable")

test_that("Binary case",{
  
  tar <- rep(0L:1, 5)
  feat <- rep(0L:1, 5)
  fast_crosstable(bit::as.bit(tar), length(tar), sum(tar),  feat)

  expect_equal(fast_crosstable(bit::as.bit(tar), length(tar), sum(tar),  feat), 
               c(5L, 0L, 0L, 5L))
})
