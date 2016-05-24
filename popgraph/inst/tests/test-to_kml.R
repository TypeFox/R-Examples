context("to_kml.R")

test_that("tests",{
  a <- matrix( 0,nrow=4,ncol=4)
  a[1,2] <- a[1,3] <- a[2,3] <- a[3,4] <-1
  a <- a + t(a)
  
  rownames(a) <- c("Olympia","St. Louis", "Ames","Richmond")
  
  graph <- as.popgraph( a )
  
  expect_that( to_kml( "bob" ), throws_error() )
  expect_that( to_kml( graph ), throws_error() )
  
  V(graph)$Latitude <- c( 47.15, 38.81, 43.08, 37.74 )
  V(graph)$Longitude <- c(-122.89,-89.98, -93.47, -77.16 )
  
  ret <- to_kml(graph)
  
  expect_that( ret, is_a("character")) 
  
  
  
})