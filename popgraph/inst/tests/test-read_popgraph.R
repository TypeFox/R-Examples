context("read.popgraph.R")

test_that("testing", {

  file <- system.file("extdata","lopho.pgraph",package="popgraph")
  if( file == "")
    file <- "/Users/rodney/Documents/Dropbox/R/popgraph/popgraph/inst/extdata/lopho.pgraph"
  
  graph <- read.popgraph( file )
  
  
  expect_that( graph, is_a("igraph") )
  expect_that( length( V(graph) ), equals(21) )
  expect_that( length( E(graph) ), equals(50) )
  expect_that( is.weighted(graph), is_true() )
  
}
)
