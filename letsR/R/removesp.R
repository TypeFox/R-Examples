# Function to remove species with values of zero (i.e. species not present in the grid) in PresenceAbscence object
# Bruno Vilela

.removeSp <- function(x) {  
  
  rem <- which(colSums(x[, -(1:2), drop = FALSE]) == 0) + 2
  
  if (length(rem) > 0) {
    x <- x[, -rem, drop = FALSE]
  }
  
  if (ncol(x) == 2) {
    stop("No species left after removing species without occurrences")
  }
  
  return(x)
}
