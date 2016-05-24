# Function to remove cells with values of zero (i.e. cells with no species present) in a PresenceAbscence object
# Bruno Vilela

.removeCells <- function(x) {
  rem <- which(rowSums(x[, -c(1, 2), drop = FALSE]) == 0)
  if (length(rem) > 0) {
    x <- x[-rem, , drop = FALSE]
  }
  if(nrow(x) == 0) {
    stop("No cells left after removing cells without occurrences")
  }
  return(x)
}
