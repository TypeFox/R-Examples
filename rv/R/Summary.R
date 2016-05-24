

Summary.rv <- function(..., na.rm=FALSE)
{
  S <- sims(c(...)) # an L x n matrix of numbers
  M <- apply(S, 1, .Generic, na.rm=na.rm)
  rvsims(M)
}
