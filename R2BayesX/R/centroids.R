centroids <- function(map) 
{
  n <- length(map)
  cp <- matrix(0, n, 2L)
  for(i in 1L:n) {
    cp[i,] <- centroidpos(na.omit(map[[i]]))
  }
  rownames(cp) <- names(map)
  colnames(cp) <- c("x", "y")

  return(cp)
}

