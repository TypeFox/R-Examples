context("Testing GitHub issue #1")


test_that("NAs are preserved in edge attributes", {
          g <- igraph::graph( c(0,1, 1,2, 2,3, 3,4, 4,2)+1, directed=TRUE)
          igraph::E(g)$label <- c(1,2,3,NA,4)
          net <- asNetwork(g)
          expect_true(any( is.na(network::get.edge.attribute(net, "label"))))
          ig <- asIgraph(net)
          expect_true(any( is.na(igraph::get.edge.attribute(ig, "label"))))
} )

test_that("NAs are preserved in vertex attributes", {
          g <- igraph::graph( c(0,1, 1,2, 2,3, 3,4, 4,2)+1, directed=TRUE)
          igraph::V(g)$label <- c(1,2,3,NA,4)
          net <- asNetwork(g)
          expect_true(any( is.na(network::get.vertex.attribute(net, "label"))))
          ig <- asIgraph(net)
          expect_true(any( is.na(igraph::get.vertex.attribute(ig, "label"))))
} )
