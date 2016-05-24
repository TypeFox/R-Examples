soft <- function(x, T) {
  y <- max(abs(x) - T, 0)
  return(y/(y+T) * x)
}

