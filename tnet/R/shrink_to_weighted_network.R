`shrink_to_weighted_network` <-
function(net){
  # Add a column
  net <- as.matrix(net)
  # Find duplicates
  net <- net[order(net[,1],net[,2]),]
  index <- !duplicated(net[,1:2])
  # Sum duplicates to get weights
  net <- cbind(net[index,], tapply(index, cumsum(index), length))
  return(as.tnet(net, type="weighted one-mode tnet"))
}