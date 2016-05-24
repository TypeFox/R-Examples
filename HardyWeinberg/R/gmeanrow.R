gmeanrow <- function(x) {
  # computes the geometric mean for a vector
  y <- exp(mean(log(x)))
  return(y)
}
