adjacency_to_shreve <- function(adjacency){
  trips        <- triplet(adjacency$adjacency)$indices
  shreve.order <- inverse.order <- vector("numeric", length = nrow(adjacency$adjacency))
  #   INVERSE ORDER PROVIDES THE NUMBER OF STREAM SEGMENTS
  #   EACH SEGMENTS LIES UPSTREAM FROM THE SOURCE
  # the root node is the segment that has no downstream neighbours, ie. the outlet
  root.node    <- which(rowSums(adjacency$adjacency) == 0)
  if(length(root.node) > 1){
    bid_roots <- nchar(adjacency$rid_bid[root.node, 2])
    root.node <- root.node[which(bid_roots == max(bid_roots))]
  }
  # inverse.order is a vector that accumulates from root.node upstream
  # each element indicates the 'height' in the network of each segment
  inverse.order[root.node] <- 1 
  # every element of inverse.order should be non-zero.  
  # A zero might indicate a lake or something disconnected from the network
  # to get the inverse.order, siply count the number of characters in the binaryID vector
  inverse.order <- nchar(adjacency$rid_bid[, 2])
  for(i in 1:nrow(trips)) inverse.order[trips[i,1]] <- inverse.order[trips[i,2]] + 1
  # sources is a vector indicating the segments that have no upstream neighbours, ie the sources of the stream
  sources <- which(colSums(adjacency$adjacency) == 0)
  # shreve.order assigns a weight of 1 to these source segments
  shreve.order[sources] <- 1
  # remaining gives the network 'heights' of the non-source segments that remain
  # in reverse order, since shreve.order will accumlate from the highest to the lowest in the network
  remaining <- rev(sort(unique(inverse.order[-sources])))
  # we now ignore the sources and simply populate the remaining shreve.orders
  inverse.order[sources] <- NA
  # need the row locations of 1's in adjacency
  # since spam uses column compression, must transpose to extract these efficiently
  tadjacency <- t(adjacency$adjacency)
  # populate elements of shreve.order in reverse order of network 'height' (larger number means 'higher')
  for(i in 1:length(remaining)){
    # identify all those segments with network height == remaining[i]
    next.down <- which(inverse.order == remaining[i])
    # cycle through each of next.down in turn and add up their upstream weights
    # this makes sense because every segment with the same network height has weight indep. of the others
    for(j in 1:length(next.down)){
      shreve.order[next.down[j]] <- sum(shreve.order[tadjacency[next.down[j],]@colindices])
    }
  }
  # finally print the shreve.order
  shreve.order
}
