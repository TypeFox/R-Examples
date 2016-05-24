
# test comparing edgelist

require(network)
require(testthat)

test<-network.initialize(5)
add.edges(test,5,1)
add.edges(test,1,5)
set.edge.attribute(test,'value',c('a','b'))
set.edge.attribute(test,'weight',10:11)

expect_equal(as.matrix.network.edgelist(test),structure(c(5L, 1L, 1L, 5L), .Dim = c(2L, 2L), n = 5, vnames = 1:5))
# sort order should be different
expect_equal(as.edgelist(test),structure(c(1L, 5L, 5L, 1L), .Dim = c(2L, 2L), n = 5, vnames = 1:5, directed = TRUE, bipartite = FALSE, loops = FALSE, inverted = FALSE, class = c("edgelist","matrix")))

expect_true(is.edgelist(as.edgelist(test)))

# numeric attribute
expect_equal(as.matrix.network.edgelist(test,attrname='weight'),structure(c(5L, 1L, 1L, 5L, 10L, 11L), .Dim = 2:3, n = 5, vnames = 1:5))

# character attribute  NOTE makes the matrix character as well
expect_equal(as.matrix.network.edgelist(test,attrname='value'),structure(c('5', '1', '1', '5', 'a', 'b'), .Dim = 2:3, n = 5, vnames = 1:5))


undir<-network.initialize(5,directed=FALSE)
add.edges(undir,5,1)
# direction will be swapped to tail < head
expect_equal(as.edgelist(undir)[,], c(1,5))

# empty network
as.edgelist(network.initialize(0))

# deleted edges
deledge<-network.initialize(5)
add.edges(deledge,1:3,2:4)
delete.edges(deledge,2)
expect_equal(as.edgelist(deledge),structure(c(1L, 3L, 2L, 4L), .Dim = c(2L, 2L), n = 5, vnames = 1:5, directed = TRUE, bipartite = FALSE, loops = FALSE, inverted = FALSE, class = c("edgelist", "matrix")))
