# This is the function to deal with a 'raw' vector
# Not exported.

# Returns NAs for non-numeric input or all-NA input

hdiVector <- function(object, credMass=0.95, ...) {
  result <- c(NA_real_, NA_real_)
  if(is.numeric(object)) {
    x <- sort(object)  # also removes NAs
    n <- length(x)
    if(n > 0) {
      # exclude <- ceiling(n * (1 - credMass)) # Not always the same as...
      exclude <- n - floor(n * credMass)       # Number of values to exclude
      low.poss <- x[1:exclude]             # Possible lower limits...
      upp.poss <- x[(n - exclude + 1):n]   # ... and corresponding upper limits
      best <- which.min(upp.poss - low.poss)      # Combination giving the narrowest interval
      result <- c(low.poss[best], upp.poss[best])
    }
  }
  names(result) <- c("lower", "upper")
  return(result)
}




