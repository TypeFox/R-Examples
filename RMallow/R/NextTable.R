#' Calculates the table of Kendall distances in (N+1)! space, given those in N!
#' space.
#' 
#' This is identical to counting the number of fully-ordered vectors at each
#' bubble sort distance in (N+1)! space.
#' 
#' 
#' @param last.table Table of distances in N! space.
#' @param N.last N
#' @return Table of distances in (N+1)! space.
#' @author Erik Gregory
#' @keywords bubblesort Kendall

NextTable <-
function(last.table, N.last) {
  len <- (N.last + 1)*(N.last)/2 + 1
  # Length of vector, minus the stuff that is easy to fill in.
  to.go <- len - 2*(N.last + 1)
  # The numbers in the next table
  nex <- c(cumsum(last.table)[1:(N.last + 1)], 
                  cumsum(as.numeric(last.table))[(N.last + 2):(N.last + 1 + to.go)] - cumsum(as.numeric(last.table))[1:to.go], 
                  rev(cumsum(as.numeric(last.table))[1:(N.last + 1)]))
  return(nex)
}
