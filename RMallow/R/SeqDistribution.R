#' Calculates distances in N! space.
#' 
#' Calculates Kendall's distances of each sequence in N! space.  This is VERY
#' Inefficient for N >= 8.  See DistanceDistribution for an astronomical
#' improvement (possibly on the order of 10^10).
#' 
#' 
#' @param N Length of the ranking.  Preferrably less than 9.
#' @return Vector of Kendall distances from 1:N to each sequence in N! space.
#' @author Erik Gregory
#' @keywords distance bubblesort

SeqDistribution <-
function(N) {
  seqs <- do.call("rbind", permn(N))
  all.dists <- AllSeqDists(seqs)
  return(all.dists)
}
