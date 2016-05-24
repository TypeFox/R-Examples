overlapTrue <-
function(d1, d2=NULL) {
#  Args:
#   d1 : either a vector or a 2-column matrix of densities
#   d2 : a vector of densities; ignored if d1 is a matrix
#   Densities should refer to equidistant points on a circle.
# Returns:
#   coefficient of overlap

  # Deal with d1 as matrix:
  if(!is.null(ncol(d1)) && ncol(d1) == 2) {
    d2 <- d1[,2]
    d1 <- d1[,1]
  }
  # Remove first value if same as the last: (added 19 Oct 2012)
  if(identical(all.equal(d1[1], d1[length(d1)]), TRUE) &&
      identical(all.equal(d2[1], d2[length(d1)]), TRUE))  {
    d1 <- d1[-1]
    d2 <- d2[-1]
  }
  # Scale each distribution to add to 1:
  d1 <- d1 / sum(d1)
  d2 <- d2 / sum(d2)
  # Compute overlap:
  sum(pmin(d1, d2))
}
