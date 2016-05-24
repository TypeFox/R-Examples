context("decorate_popgraph.R")

test_that("testing", {
  
  A <- matrix(0, nrow=4, ncol=4)
  A[1,2] <- A[2,3] <- A[1,3] <- A[3,4] <- 1
  A <- A + t(A)
  rownames(A) <- colnames(A) <- LETTERS[1:4]
  g <- as.popgraph( A )

  df <- data.frame(name=LETTERS[1:4], Position=c(1,2,3,4))  
  
  expect_that(decorate_graph( g, df ), throws_error() )
  
  g.df <- decorate_graph( g, df, stratum="name" )
      
  expect_that( g.df, is_a("igraph") )
  expect_that( g.df, is_a("popgraph"))

  expect_that( V(g.df)$Bob, is_a("NULL") )
  
  expect_that( V(g.df)$Position, is_a("numeric"))  
  expect_that( V(g.df)$Position[1], equals(1) )
}
)
