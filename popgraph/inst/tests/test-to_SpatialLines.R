context("to_SpatialLines.R")

test_that("testing", {

  a <- matrix( 0,nrow=4,ncol=4)
  a[1,2] <- a[1,3] <- a[2,3] <- a[3,4] <-1
  a <- a + t(a)
  
  rownames(a) <- c("Olympia","St. Louis", "Ames","Richmond")
  
  graph <- as.popgraph( a )
  
  expect_that( to_SpatialLines( "bob" ), throws_error() )
  expect_that( to_SpatialLines( graph ), throws_error() )

               
  V(graph)$Latitude <- c( 47.15, 38.81, 43.08, 37.74 )
  V(graph)$Longitude <- c(-122.89,-89.98, -93.47, -77.16 )
               
  sl <- to_SpatialLines( graph )
  
  expect_that( sl, is_a("SpatialLines") )
  bbexp <- matrix(c(-122.89,-77.16,37.74,47.15),ncol=2, byrow=T)
  expect_that( bbox(sl), is_equivalent_to(bbexp))
  
}
)
