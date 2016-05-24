context("overlay_popgraph.R")

test_that("tests", {
  
  A <- matrix(0, nrow=5, ncol=5)
  A[1,2] <- A[2,3] <- A[1,3] <- A[3,4] <- A[4,5] <- 1
  A <- A + t(A)
  
  g <- as.popgraph( A )
  
  #require(maps)
  #map("state")
  #overlay_popgraph( g )  
  
  expect_that( overlay_popgraph(FALSE), throws_error() )
  V(g)$Longitude <- c(-122.89,-122.49,-89.98, -93.47, -77.16 )
  expect_that( overlay_popgraph(graph), throws_error())
  V(g)$Latitude <- c( 47.15, 48.75,38.81, 42.26, 37.74 )
  expect_that( overlay_popgraph(graph), throws_error())  
  V(g)$name <- c("Olympia","Bellingham","St. Louis","Ames","Richmond")
  expect_that( overlay_popgraph(graph), throws_error())
  
})