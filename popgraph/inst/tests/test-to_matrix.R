context("to_matrix.R")

test_that("testing", {
  
  a <- matrix( 0,nrow=4,ncol=4)
  a[1,2] <- a[1,3] <- a[2,3] <- a[1,4] <-1
  a <- a + t(a)
  g <- as.popgraph( a )
  
  K <- length(V(g))
  m1 <- to_matrix( g, mode="adjacency")
  m2 <- to_matrix( g, mode="shortest path")
  m3 <- to_matrix( g, mode="edge weight")
  
  # type
  expect_that(m1, is_a("matrix"))
  expect_that(m2, is_a("matrix"))
  expect_that(m3, is_a("matrix"))
  
  # size
  expect_that(dim(m1), is_equivalent_to(c(K,K)))
  expect_that(dim(m2), is_equivalent_to(c(K,K)))
  expect_that(dim(m3), is_equivalent_to(c(K,K)))
  
  #diagonal
  expect_that(diag(m1), is_equivalent_to(rep(0,K)))
  expect_that(diag(m2), is_equivalent_to(rep(0,K)))
  expect_that(diag(m3), is_equivalent_to(rep(0,K)))
}
)