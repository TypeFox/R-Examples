context("Testing whether custom 'attrmap' specifications work")

test_that("'na' vertex attribute can be dropped when going network->igraph", {
          net <- asNetwork(exIgraph)
          a <- attrmap()
          r <- data.frame(type="vertex", fromcls="network",
                          fromattr="na", tocls="igraph", toattr=NA,
                          stringsAsFactors=FALSE)
          aa <- rbind(a, r)
          g <- asIgraph(net, amap=aa)
          expect_false( "na" %in% igraph::list.vertex.attributes(g))
} )
