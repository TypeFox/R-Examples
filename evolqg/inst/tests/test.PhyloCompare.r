test_that("PhyloCompare gives correct results", 
{
  TREE <- "iris((setosa:1,versicolor:1):1,virginica:2);"
  tree.iris <- ape::read.tree(text = TREE)
  iris.cov.list <- dlply(iris, "Species", function(x) cov(x[,1:4]))
  cov.matrices <- PhyloW(tree.iris, iris.cov.list)
  RS.iris = PhyloCompare(tree.iris, cov.matrices)
  mantel.iris = PhyloCompare(tree.iris, llply(cov.matrices, cov2cor), function(x, y) cor(x[lower.tri(x)],
                                                                                         y[lower.tri(y)]))
  krz.iris = PhyloCompare(tree.iris, cov.matrices, KrzCor)
  dist.iris = PhyloCompare(tree.iris, cov.matrices, MatrixDistance, "Over")
  
  expect_that(RS.iris, is_a('list'))
  expect_that(mantel.iris, is_a('list'))
  expect_that(krz.iris, is_a('list'))
  expect_that(dist.iris, is_a('list'))
  
  expect_that(RS.iris[[1]], is_a('data.frame'))
  expect_that(mantel.iris[[1]], is_a('data.frame'))
  expect_that(krz.iris[[1]], is_a('data.frame'))
  expect_that(dist.iris[[1]], is_a('data.frame'))
  
  expect_that(RS.iris[[2]], is_a('list'))
  expect_that(mantel.iris[[2]], is_a('list'))
  expect_that(krz.iris[[2]], is_a('list'))
  expect_that(dist.iris[[2]], is_a('list'))
  
  internal.check = RS.iris[[1]][,2] == llply(unique(RS.iris[[1]]$node), function(node) unique(unlist(RS.iris[[2]][(grep(paste0(node, ','), names(RS.iris[[2]])))])))
  expect_true(all(internal.check))
  
  internal.check = mantel.iris[[1]][,2] == llply(unique(mantel.iris[[1]]$node), function(node) unique(unlist(mantel.iris[[2]][(grep(paste0(node, ','), names(mantel.iris[[2]])))])))
  expect_true(all(internal.check))
  
  internal.check = krz.iris[[1]][,2] == llply(unique(krz.iris[[1]]$node), function(node) unique(unlist(krz.iris[[2]][(grep(paste0(node, ','), names(krz.iris[[2]])))])))
  expect_true(all(internal.check))
  
  internal.check = dist.iris[[1]][,2] == llply(unique(dist.iris[[1]]$node), function(node) unique(unlist(dist.iris[[2]][(grep(paste0(node, ','), names(dist.iris[[2]])))])))
  expect_true(all(internal.check))
})