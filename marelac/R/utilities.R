
## This function checks if the length of all vectors in a list
## is either "1" or has the same maximum value
checkVecLength <- function(x) {
  if (!is.list(x)) stop("x must be a list of vectors")
  alleng <- sapply(x, length)
  maxlen <- max(alleng)
  all((alleng == maxlen) | (alleng == 1))
}
