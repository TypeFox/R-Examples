# (c) Kevin Dunn, 2015.

popcorn <- function(T=120){
  # Simulates a stovetop popcorn cooking process where there is 1 factor only:
  # cooking time, T, measured in seconds.
  # The same number of kernels are cooked each time using the same heat setting on the stove.
  # 
  # The outcome is: number of popped kernels from the batch.
  
  if (length(T) > 1){
    stop("Cooking popcorn batches in parallel is (intentionally) not allowed.")
  }
  if (!is.finite(T)){
    stop("Please provide finite numeric values as inputs.")
  }
  if (T < 77) {
    stop("No popcorn was made: please cook for longer.")
  } else {
    coded <- (T - 135.0) / 15.0
    y <- coded * 15 - 2.4 * coded * coded + 93 + (runif(1, min = 0, max = 1) * 6 - 3.0)
    y <- max(0, round(y))
  }
  return(y)
}