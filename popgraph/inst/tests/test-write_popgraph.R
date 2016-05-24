context("write_popgraph.R")

test_that("testing", {

  expect_that( write_popgraph(FALSE), throws_error() )

  a <- matrix( 0,nrow=4,ncol=4)
  a[1,2] <- a[1,3] <- a[2,3] <- a[1,4] <-1
  a <- a + t(a)
  graph <- as.popgraph(a)
  
  expect_that( write_popgraph(graph), throws_error() )
}
)
