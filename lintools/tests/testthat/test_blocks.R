
context("blocks")


test_that("block index",{
  A <- matrix(c(
    1,1,0,
    0,1,0,
    0,0,1),byrow=TRUE,nrow=3)

  expect_equal(block_index(A)
      ,list(c(1,2),3)
  )
})


