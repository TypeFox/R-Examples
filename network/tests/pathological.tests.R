library(network)
if (require(statnet.common,quietly=TRUE)){

  opttest({
    gctorture(TRUE)
    n <- 10
    test <- network.initialize(n)
    for (i in 1:(n-1)){
      for (j in (i+1):n){
        cat(i,j,'\n')
        get.inducedSubgraph(test,v=i:j)
      }
    }
    gctorture(FALSE)
  },'Ticket #180 Test 1','NETWORK_pathology_TESTS')

  opttest({
    gctorture(TRUE)
    test <- network.initialize(10)
    delete.vertices(test,5)
    gctorture(FALSE)
  },'Ticket #180 Test 2','NETWORK_pathology_TESTS')

  opttest({
    x <- network.initialize(10)
    x[,] <- 1
    try(set.edge.value(x,'foo',matrix('bar',5,5)))
  },'Ticket #827','NETWORK_pathology_TESTS')

}
