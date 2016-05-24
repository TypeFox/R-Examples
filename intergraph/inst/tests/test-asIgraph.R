context("Creating igraphs from data frames")


test_that("Disassembling to d.f and assembling back to igraph gives the same result", {
  # convert to data frames
  l <- asDF(exIgraph)
  # assemble back
  g <- asIgraph( l$edges, vertices=l$vertexes)
  expect_true( igraph::identical_graphs(g, exIgraph) )
} )



test_that("Providing non-existing name to 'vnames' throws an error", {
  ## testing 'vnames' argument
  # non-existent column in 'vertices'
  expect_that( asIgraph(l$edges, vertices=l$vertexes, vnames="foo"), throws_error())
} )



test_that("Vertex names are properly set via 'vnames' argument for directed network", {
  # existing column in 'vertices'
  l <- asDF(exIgraph)
  g <- asIgraph( l$edges, vertices=l$vertexes, vnames="label")
  expect_equal( l$vertexes$label, igraph::get.vertex.attribute(g, "name"))
} )



test_that("Vertex names are properly set via 'vnames' argument for undirected network", {
  ## above tests but for the undirected network
  ## convert to data frames and assemble back to igraph object
  l <- asDF(exIgraph2)
  g <- asIgraph( l$edges, vertices=l$vertexes, directed=FALSE)
  expect_true( igraph::identical_graphs(g, exIgraph2) )
} )










context("Creating igraphs from networks")



test_that("Conversion for exNetwork is OK tested with netcompare", {
  # directed network
  res <- netcompare( asIgraph(exNetwork), exNetwork, test=TRUE )
  expect_true(res)
} )

test_that("Conversion for exNetwork2 is OK tested with netcompare", {
  # undirected network
  res2 <- netcompare( asIgraph(exNetwork2), exNetwork2, test=TRUE )
  expect_true(res2)
} )


test_that("Conversion of bipartite networks is not yet supported", {
  ### bipartite network (not yet supported)
  m <- matrix(0, ncol=2, nrow=3)
  m[1,1] <- m[2,1] <- m[1,2] <- m[3,2] <- 1
  net <- network::network(t(m), bipartite=TRUE, directed=FALSE)
  expect_error(asIgraph(net))
} )
