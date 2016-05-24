context("congruence_topology.R")

test_that("testing", {

  A <- matrix(0,nrow=4,ncol=4)
  B <- A
  
  A[1,2] <- A[1,3] <- A[2,3] <- A[3,4] <- 1
  B[1,2] <- B[1,3] <- B[1,4] <- B[2,4] <- B[2,3] <- 1
  
  A <- A + t(A)
  B <- B + t(B)
  
  rownames(A) <- colnames(A) <- rownames(B) <- colnames(B) <- LETTERS[1:4]
  
  graph1 <- as.popgraph(A)
  graph2 <- as.popgraph(B)
  
  cong <- congruence_topology(graph1,graph2)
  expect_that( cong, is_a("igraph") )
  expect_that( cong, is_a("popgraph"))
  expect_that( length(V(cong)), equals(4) )
  expect_that( length(E(cong)), equals(3) )
}
)
