re_map_rid <- function(rid_vector, all_rid){
  # in cases where multiple networks are present, gaps could occur in the rid sequence
  # this function maps the rid vector in a simple way so that it agrees with the adjacency matrix
  mapped_rid <- vector("numeric", length = length(rid_vector))
  rng_rid    <- range(all_rid)
  new_key <- 1:length(unique(all_rid))
  old_rid <- sort(all_rid)
  for(i in 1:length(rid_vector)){
    mapped_rid[i] <- new_key[which(old_rid == rid_vector[i])]
  }
  mapped_rid
}