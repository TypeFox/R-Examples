context("to_pgraph.R")

test_that("tests", {
  a <- matrix( 0,nrow=4,ncol=4)
  a[1,2] <- a[1,3] <- a[2,3] <- a[1,4] <- 2
  a <- a + t(a)
  rownames(a) <- colnames(a) <- as.character(LETTERS[1:4])
  
  graph <- as.popgraph( a )
  
  expect_that( to_pgraph(FALSE), throws_error() )

  ret <- to_pgraph( graph )
  expect_that( ret, is_equivalent_to("4\t4\nA\t1\t1\nB\t1\t1\nC\t1\t1\nD\t1\t1\nA\tB\t2\nA\tC\t2\nA\tD\t2\nB\tC\t2"))
})
