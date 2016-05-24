
check_additive_weight <- function(adjacency, weight){
  # transpose the adjacency so that the rows are accessible (as cols in the transpose)
  tadjacency <- t(adjacency)
  # create binary numeric vector, i^th element is a check that weights immediately 
  # upstream to i sum exactly to the weight of i
  is.additive <- vector("numeric", length = nrow(adjacency))
  # create a numeric vector that measures the discrepancy between upstream and downstream weights
  discrep <- vector("numeric", length = nrow(adjacency))
  # loop through the rows of the adjacency matrix, performing the sum  check for each
  for(i in 1:nrow(adjacency)){
    dn.weight <- weight[i]
    col.inds  <- tadjacency[i,]
    if(!sum(col.inds) == 0){
      up.weight <- sum(weight[col.inds@colindices])
      is.additive[i] <- as.numeric(!dn.weight == up.weight)
      discrep[i] <- sum(dn.weight - up.weight)^2
    } else {
      is.additive[i] <- 0
    }
  }
  # it is possible that the additivity condition holds but is subject to round off and other small errors
  # therefore only print a warning if the rmse reaches a sensible tolerance
  if((mean(is.additive)) > 0 & (mean(discrep) > 10^-5)){
    rmse <- mean(discrep)
    commt <- paste("Additivity of selected weight variable doesn't seem to hold, \nRMSE = ", rmse, 
                   " was observed. \nProceed with caution.", sep = "")
    warning(commt)
  }
}