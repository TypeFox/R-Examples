library(bnstruct)
context("Testing CPTs")

dataset <- child()
dataset <- impute(dataset)
net <- learn.network(dataset)

dag <- dag(net)
s <- sapply(1:num.nodes(net), function(x) {
  parents <- which(dag[,x] > 0)
  prod(node.sizes(net)[parents])
})

test_that("sum of cpts", {
  expect_equal(sum(unlist(cpts(net))), sum(s))
})

