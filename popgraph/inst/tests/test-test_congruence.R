context("test_congruence.R")

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
  
  ret <- test_congruence(graph1,graph2)
  expect_that( ret, is_a("htest"))
  expect_that( ret$parameter, is_equivalent_to(4) )
  
  
  ret <- test_congruence(graph1,graph2,method="combinatorial")
  expect_that( ret, is_equivalent_to(2/3))
  
  
  expect_that( test_congruence(graph1,graph2,method="structural"), throws_error())
  
  
  A <- matrix(0,nrow=4,ncol=4)
  graphA <- graph.adjacency(A,mode="undirected")
  expect_that( test_congruence(graph1,graphA), throws_error())
  
  
  
}
)
