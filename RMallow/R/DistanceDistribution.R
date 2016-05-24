#' Calculate the Kendall distance distribution in N! space.
#' 
#' This function counts the number of fully-ordered vectors at each distance in
#' N! space.
#' 
#' 
#' @param N Integer value, greater than or equal to 3.
#' @return Table-like structure, where the names represent the distance from
#' the modal sequence of each sequence in N! space, and the values represent the number of sequences at that distance in
#' the sequence space.
#' @author Erik Gregory
#' @keywords bubblesort Kendall
#' # DistanceDistribution(10)

DistanceDistribution <-
function(N = 3) {
  if (N <= 4) {
  	out <- table(SeqDistribution(N))
  }
  else {
  	out <- table(SeqDistribution(4))
        i <- 4
	while (i < N) {
	   out <- NextTable(out, i)
	   i <- i + 1
	}
  	names(out) <- as.character(0:(length(out) - 1))
  }
  return(out)
}
