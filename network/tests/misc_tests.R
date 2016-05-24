# tests for misc R functions

require(network)
require(testthat)
# tests for has.edges

test<-network.initialize(5)
test[1,2]<-1
expect_equal(has.edges(test), c(TRUE,TRUE,FALSE,FALSE,FALSE))
expect_equal(has.edges(test,v=2:3),c(TRUE,FALSE))
expect_error(has.edges(test,v=10),regexp = 'argument must be a valid vertex id')
expect_equal(length(has.edges(network.initialize(0))),0)
