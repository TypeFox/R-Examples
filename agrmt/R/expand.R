expand <- function(F) {
  # Converts frequency vector in a population vector
  # Argument: F = frequency vector
  # Example: F <- c(10,20,34,5)
  # input validation: is this a frequency vector?
  if (min(F) < 0) stop("Error: negative values found in frequency vector.")
  if (all(round(F,0) == F) == FALSE) stop("Error: not all values are whole numbers; not a frequency vector; multiply by 100?")

  k <- length(F) # prepare expansion
  D <- NULL      # prepare expansion
  for (i in 1:k) {
    D <- c(D,rep(i,F[i]))   # expand
  }
  return(D)      # population vector
}
