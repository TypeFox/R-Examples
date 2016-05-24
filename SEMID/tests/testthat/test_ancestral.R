library(SEMID)
context("Components related to generic identifiability by ancestor decomposition.")

rUndirectedAdjMat = function(n, p) {
  mat = matrix(runif(n^2), ncol=n) < p
  mat = mat * upper.tri(mat)
  return(mat + t(mat))
}
rConnectedAdjMatrix = function(n, p) {
  weights = runif(n*(n-1)/2)
  g = minimum.spanning.tree(graph.full(n), weights=weights)
  adjMatrix = as.matrix(get.adjacency(g))
  adjMatrix = (upper.tri(matrix(0,n,n)) & matrix(sample(c(T, F), n^2, replace=T, prob=c(p, 1-p)), ncol=n)) | adjMatrix
  adjMatrix = 1*(adjMatrix | t(adjMatrix))
  return(adjMatrix)
}
rDirectedAdjMatrix = function(n,p) {
  return(1*(upper.tri(matrix(0,n,n)) & matrix(sample(c(T, F), n^2, replace=T, prob=c(p, 1-p)), ncol=n)))
}
rDirectedAdjMatrix = function(n, p) {
  return(1*(upper.tri(matrix(0,n,n)) & matrix(sample(c(T, F), n^2, replace=T, prob=c(p, 1-p)), ncol=n)))
}
getAdjMat = function(g) { as.matrix(get.adjacency(g)) }

test_that("ancestors function works as expected.", {
  ## Output should be empty
  expect_equal(ancestors(graph.empty(), c()), numeric(0))
  expect_equal(ancestors(graph.empty(), c(1,2,3)), numeric(0))
  expect_equal(ancestors(graph.star(5, mode="in"), c()), numeric(0))

  # Graph with no edges
  expect_equal(ancestors(graph.empty(5), 1), 1)
  expect_equal(ancestors(graph.empty(5), 2), 2)
  expect_equal(ancestors(graph.empty(5), c(2,4,5)), c(2,4,5))

  # Graphs with edges
  expect_equal(ancestors(graph.star(5, mode="in"), 1), 1:5)
  expect_equal(ancestors(graph.star(5, mode="in"), 2), 2)
  expect_equal(ancestors(graph.star(5, mode="out"), 1), 1)
  expect_equal(ancestors(graph.full(5, directed=T), 1), 1:5)
  g = graph.edgelist(matrix(c(1,2,2,3,3,5,5,7,2,4,4,6,6,7,4,5), ncol=2, byrow=T))
  expect_equal(ancestors(g, 1), 1)
  expect_equal(ancestors(g, 2), c(1,2))
  expect_equal(ancestors(g, 4), c(1,2,4))
  expect_equal(ancestors(g, 5), c(1,2,3,4,5))
  expect_equal(ancestors(g, 6), c(1,2,4,6))
  expect_equal(ancestors(g, c(6, 3)), c(1,2,3,4,6))
  expect_equal(ancestors(g, 7), 1:7)
})

test_that("parents function works as expected.", {
  ## Output should be empty
  expect_equal(parents(graph.empty(), c()), numeric(0))
  expect_equal(parents(graph.empty(), c(1,2,3)), numeric(0))
  expect_equal(parents(graph.star(5, mode="in"), c()), numeric(0))

  # Graph with no edges
  expect_equal(parents(graph.empty(5), 1), 1)
  expect_equal(parents(graph.empty(5), 2), 2)
  expect_equal(parents(graph.empty(5), c(2,4,5)), c(2,4,5))

  # Graphs with edges
  expect_equal(parents(graph.star(5, mode="in"), 1), 1:5)
  expect_equal(parents(graph.star(5, mode="in"), 2), 2)
  expect_equal(parents(graph.star(5, mode="out"), 1), 1)
  expect_equal(parents(graph.full(5, directed=T), 1), 1:5)
  g = graph.edgelist(matrix(c(1,2,2,3,3,5,5,7,2,4,4,6,6,7,4,5,1,7), ncol=2, byrow=T))
  expect_equal(parents(g, 1), 1)
  expect_equal(parents(g, 2), c(1,2))
  expect_equal(parents(g, 4), c(2,4))
  expect_equal(parents(g, 5), c(3,4,5))
  expect_equal(parents(g, 6), c(4,6))
  expect_equal(parents(g, c(6, 3)), c(2,3,4,6))
  expect_equal(parents(g, 7), c(1,5,6,7))
})

test_that("siblings function works as expected.", {
  ## Output should be empty
  expect_equal(siblings(graph.empty(), c()), numeric(0))
  expect_equal(siblings(graph.empty(), c(1,2,3)), numeric(0))
  expect_equal(siblings(graph.star(5, mode="in"), c()), numeric(0))

  # Graph with no edges
  expect_equal(siblings(graph.empty(5), 1), 1)
  expect_equal(siblings(graph.empty(5), 2), 2)
  expect_equal(siblings(graph.empty(5), c(2,4,5)), c(2,4,5))

  # Graphs with edges
  expect_equal(siblings(graph.star(5, mode="in"), 1), 1:5)
  expect_equal(siblings(graph.star(5, mode="in"), 2), c(1,2))
  expect_equal(siblings(graph.star(5, mode="out"), 1), 1:5)
  expect_equal(siblings(graph.star(5, mode = "undirected"), 1), 1:5)
  expect_equal(siblings(graph.full(5, directed=F), 1), 1:5)
  g = graph.edgelist(matrix(c(1,2,2,3,3,5,5,7,2,4,4,6,6,7,4,5,1,7), ncol=2, byrow=T), directed = F)
  expect_equal(siblings(g, 1), c(1,2,7))
  expect_equal(siblings(g, 2), c(1,2,3,4))
  expect_equal(siblings(g, c(6, 3)), c(2,3,4,5,6,7))
})

test_that("getMixedCompForNode function works as expected.", {
  ## Graph with single node
  dG = graph.empty(1, directed=T)
  bG = graph.empty(1, directed=F)
  expect_error(getMixedCompForNode(dG, bG, 1, 1)) # Since vertices are unnamed
  V(dG)$names = 1
  V(bG)$names = 1

  compList = getMixedCompForNode(dG, bG, 1, 1)
  expect_equal(compList, list(biNodes=1, inNodes=numeric(0)))
  expect_error(getMixedCompForNode(dG, bG, 1, 2))
  expect_error(getMixedCompForNode(dG, bG, 2, 1))

  # Graph with 2 nodes
  dG = graph.empty(2, directed=T)
  bG = graph.empty(2, directed=F)
  V(dG)$names = 1:2
  V(bG)$names = 2:1
  expect_error(getMixedCompForNode(dG, bG, 1, 1)) # Since vertices are wrongly named

  V(dG)$names = 1:2
  V(bG)$names = 1:2

  compList = getMixedCompForNode(dG, bG, 1, 1)
  expect_equal(compList, list(biNodes=1, inNodes=numeric(0)))

  compList = getMixedCompForNode(dG, bG, 1, 2)
  expect_equal(compList, list(biNodes=numeric(0), inNodes=numeric(0)))

  compList = getMixedCompForNode(dG, bG, 2, 1)
  expect_equal(compList, list(biNodes=numeric(0), inNodes=numeric(0)))

  # The sink marginalization example
  dG = graph.edgelist(matrix(c(1,2, 1,3, 1,6, 2,3, 2,4, 2,5, 2,6, 3,4, 4,5),
                             ncol=2, byrow=T))
  bG = graph.edgelist(matrix(c(1,6, 1,4, 2,3, 2,5, 2,6),
                             ncol=2, byrow=T), directed=F)
  V(dG)$names = 1:6
  V(bG)$names = 1:6
  compList = getMixedCompForNode(dG, bG, 1:6, 1)
  expect_equal(compList, list(biNodes=1:6, inNodes=numeric(0)))
  compList = getMixedCompForNode(dG, bG, c(4,5), 5)
  expect_equal(compList, list(biNodes=5, inNodes=4))
  compList = getMixedCompForNode(dG, bG, 1:5, 5)
  expect_equal(compList, list(biNodes=c(2,3,5), inNodes=c(1,4)))
  compList = getMixedCompForNode(dG, bG, c(2,3,4), 4)
  expect_equal(compList, list(biNodes=c(4), inNodes=c(2,3)))
  compList = getMixedCompForNode(dG, bG, rev(c(2,4,5,6)), 2)
  expect_equal(compList, list(biNodes=c(2,5,6), inNodes=c(4)))

  bG = delete.edges(bG, get.edge.ids(bG, c(1,6)))
  compList = getMixedCompForNode(dG, bG, 1:6, 1)
  expect_equal(compList, list(biNodes=c(1,4), inNodes=c(2,3)))
  compList = getMixedCompForNode(dG, bG, 1:6, 3)
  expect_equal(compList, list(biNodes=c(2,3,5,6), inNodes=c(1,4)))
})

test_that("getMaxFlow function works as expected.", {
  ## Graph with single node
  dG = graph.empty(1, directed=T)
  bG = graph.empty(1, directed=F)
  expect_equal(getMaxFlow(getAdjMat(dG), getAdjMat(bG), c(), c(1), numeric(0), 1), 0)

  ## Graph with two disconnected nodes
  dG = graph.empty(2, directed=T)
  bG = graph.empty(2, directed=F)
  expect_error(getMaxFlow(getAdjMat(dG), getAdjMat(bG), c(1), c(1), numeric(0), 1))
  expect_error(getMaxFlow(getAdjMat(dG), getAdjMat(bG), c(2), c(2), numeric(0), 2))
  expect_equal(getMaxFlow(getAdjMat(dG), getAdjMat(bG), c(2), c(1), numeric(0), 1), 0)
  expect_equal(getMaxFlow(getAdjMat(dG), getAdjMat(bG), c(), c(2), numeric(0), 2), 0)
  expect_equal(getMaxFlow(getAdjMat(dG), getAdjMat(bG), c(), c(1), numeric(0), 1), 0)
  expect_error(getMaxFlow(getAdjMat(dG), getAdjMat(bG), c(), c(1), numeric(0), 2))

  ## Two connected nodes
  dG = graph.edgelist(t(c(1,2)), directed=T)
  bG = graph.empty(2, directed=F)
  expect_equal(getMaxFlow(getAdjMat(dG), getAdjMat(bG), c(), c(2), numeric(0), 2), 0)
  expect_equal(getMaxFlow(getAdjMat(dG), getAdjMat(bG), c(1), c(2), numeric(0), 2), 0)
  expect_equal(getMaxFlow(getAdjMat(dG), getAdjMat(bG), c(1), c(2), c(1), 2), 1)

  bG = add.edges(bG, c(1,2))
  expect_error(getMaxFlow(getAdjMat(dG), getAdjMat(bG), c(), c(2), c(1), 2))
  expect_error(getMaxFlow(getAdjMat(dG), getAdjMat(bG), c(1), c(1,2), c(), 2))
  expect_error(getMaxFlow(getAdjMat(dG), getAdjMat(bG), c(2), c(2), c(), 2))
  expect_equal(getMaxFlow(getAdjMat(dG), getAdjMat(bG), c(), c(2), c(), 2), 0)

  dG = add.edges(add.vertices(dG, 1), c(3,1))
  bG = add.vertices(bG, 1)
  expect_equal(getMaxFlow(getAdjMat(dG), getAdjMat(bG), c(3), c(1,2), c(3), 2), 1)
  expect_equal(getMaxFlow(getAdjMat(dG), getAdjMat(bG), c(3), c(1,2), c(3), 1), 1)
  expect_equal(getMaxFlow(getAdjMat(dG), getAdjMat(bG), c(3), c(1), c(3), 1), 1)
  expect_error(getMaxFlow(getAdjMat(dG), getAdjMat(bG), c(3), c(2,3), c(1), 2))

  ## The sink marginalization example
  dG = graph.edgelist(matrix(c(1,2, 1,3, 1,6, 2,3, 2,4, 2,5, 2,6, 3,4, 4,5),
                             ncol=2, byrow=T))
  bG = graph.edgelist(matrix(c(1,6, 1,4, 2,3, 2,5, 2,6),
                             ncol=2, byrow=T), directed=F)
  L = getAdjMat(dG)
  O = getAdjMat(bG)
  expect_equal(getMaxFlow(L, O, allowedNodes=c(1,3), biNodes=1:6, inNodes=c(), node=5), 2)
  expect_equal(getMaxFlow(L, O, allowedNodes=c(1), biNodes=1:6, inNodes=c(), node=5), 1)
  expect_equal(getMaxFlow(L, O, allowedNodes=c(6), biNodes=1:6, inNodes=c(), node=5), 1)
  L1 = L; L1[1,3] = 0
  O1 = O; O1[2,c(3,5)] = 0; O1[c(3,5), 2] = 0
  expect_equal(getMaxFlow(L1, O1, allowedNodes=c(2,6), biNodes=1:6, inNodes=c(), node=5), 1)
})

test_that("graphID.ancestralID function works as expected.", {
  # Random test
  set.seed(23231)
  ps = c(.2, .4, .6, .8)
  sims = 10
  ns = c(2, 4, 6)
  for(p in ps) {
    for(n in ns) {
      for(i in 1:sims) {
        L = rDirectedAdjMatrix(n, p)
        O = rConnectedAdjMatrix(n, p)
        gia = graphID.ancestralID(L, O)
        gig = graphID.htcID(L, O)
        gin = graphID.nonHtcID(L, O)
        expect_true(length(gia) >= length(gig))
        expect_true(all(gig %in% gia))
        expect_true((length(gia) != n) || !gin)
      }
    }
  }

  dG = graph.edgelist(matrix(c(1,2, 1,3, 1,6, 2,3, 2,4, 2,5, 2,6, 3,4, 4,5),
                             ncol=2, byrow=T))
  bG = graph.edgelist(matrix(c(1,6, 1,4, 2,3, 2,5, 2,6),
                             ncol=2, byrow=T), directed=F)
  expect_equal(as.numeric(sort(graphID.ancestralID(getAdjMat(dG), getAdjMat(bG)))), 1:6)
})
