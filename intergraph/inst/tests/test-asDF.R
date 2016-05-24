context("Converting directed networks to data frames")

test_that("Created edgelist has appropriate number of edges", {
	l <- asDF(exNetwork)
  expect_equal( nrow(l$edges), network::network.edgecount(exNetwork))
} )

test_that("Created vertex data frame has correct number of vertices", {
	l <- asDF(exNetwork)
  expect_equal( nrow(l$vertexes), network::network.size(exNetwork))
} )






context("Converting undirected networks to data frames")

test_that("Created edgelist has appropriate number of edges", {
	l <- asDF(exNetwork2)
  expect_equal( nrow(l$edges), network::network.edgecount(exNetwork))
} )

test_that("Created vertex data frame has correct number of vertices", {
	l <- asDF(exNetwork2)
  expect_equal( nrow(l$vertexes), network::network.size(exNetwork))
} )
