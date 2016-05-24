`rand.index` <- 
function(group1, group2) {
  #which points are in the same vs different groups? different groups == 1, same == 0
  x <- abs(sapply(group1, function(x) x-group1))
  x[x > 1] <- 1
  y <- abs(sapply(group2, function(x) x-group2))
  y[y > 1] <- 1
  #what pairs share the same relationships between group1 and group2 (same/same or different/different)
  sg <- sum(abs(x-y))/2
  bc <- choose(dim(x)[1],2)
  ri <- 1-sg/bc
  return(ri)
}
