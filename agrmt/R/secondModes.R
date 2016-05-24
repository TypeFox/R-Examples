secondModes <-
function(V, pos=FALSE, tolerance=0.1) {
  m1 <- modes(V, tolerance=tolerance, pos=pos)
  V2 <- V[-m1$positions]
  m2 <- modes(V2, tolerance=tolerance, pos=pos)
  a <- list(m1$at, m2$at)
  f <- list(m1$frequencies, m2$frequencies)
  m <- list(m1$mode, m2$mode)
  p <- list(m1$positions, m2$positions)
  c <- list(m1$contiguous, m2$contiguous)
  r <- list(at = a, frequencies = f, mode = m, positions = p, contiguous = c)
  # compare output to modes()
  return(r)
  }
