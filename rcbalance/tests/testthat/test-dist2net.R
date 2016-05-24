library(rcbalance)
library(optmatch)
context('Converting distance to network')

data(nuclearplants)
#reorder nuclearplants dataframe so treated units come first
nuclearplants <- nuclearplants[order(nuclearplants$pr, decreasing = TRUE),]

#switch positions of last couple of controls
nuclearplants[31:32,] <- nuclearplants[32:31,]

sparse.dist <- build.dist.struct(z = nuclearplants$pr, 
	X = subset(nuclearplants[c('date','t1','t2','cap','bw','ne','cum.n')]), 
	exact = nuclearplants$cap)

test_that('Controls with no partners are not dropped',{
	my.net <- dist2net(sparse.dist, k = 1, exclude.treated = TRUE, ncontrol = sum(nuclearplants$pr == 0))
	expect_equal(length(my.net$b) - 1 , nrow(nuclearplants))
})

test_that('Compatible with InfinitySparseMatrix and BlockedInfinitySparseMatrix', {
  sparse.mat <- match_on(pr ~ date  + t1 + t2  + cap + bw + ne, data = nuclearplants, caliper = 4)
  block.mat <- match_on(pr ~ date  + t1 + t2  + cap + bw + ne, data = nuclearplants) + exactMatch(pr ~ pt, data = nuclearplants)
  sparse.net <- dist2net.matrix(sparse.mat, k = 1)
  expect_equal(length(sparse.net$b)-1, sum(dim(sparse.mat)))
  expect_equal(sparse.net$tcarcs, sum(is.finite(sparse.mat)))
  block.net <- dist2net.matrix(block.mat, k =1)
  expect_equal(length(block.net$b)-1, sum(dim(block.mat)))
  expect_equal(block.net$tcarcs, sum(is.finite(block.mat)))  	
})


#TODO: make sure I don't round distances too sloppily.  Compare to optmatch
#TODO: performance benchmarking for large problems with InfinitySparseMatrix and BlockedInfinitySparseMatrix

#TODO (not in this file though): fix up message when you load rcbalance to make dependence on optmatch more explicit.