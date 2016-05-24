

range.rv <- function(..., na.rm=FALSE, finite=FALSE) {
  sims <- sims(c(...)) # an L x n matrix of numbers
  m <- apply(sims, 1, 'range', na.rm=na.rm, finite=finite) # gives a 2xL matrix!!!
  r <- rvsims(t(m)) # Transpose to get an L x 2 matrix.
  names(r) <- c('min', 'max')
  return(r)
}
