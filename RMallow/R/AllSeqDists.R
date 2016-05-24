#' Calculate all distances between a set of sequences and a fixed sequence.
#' 
#' Used to calculate the sequence Kendall distance distribution in N! space.
#' 
#' 
#' @param seqs Matrix or data frame of sequences.
#' @return Vector of the distances from the sequences to 1:N.
#' @author Erik Gregory
#' @keywords Kendall Distance


AllSeqDists <-
function(seqs) {
  modal <- 1:ncol(seqs)
  infos <- KendallInfo(seqs)
  dists <- apply(infos, 1, function(i) length(which(i == 1)))
  return(dists)
}
