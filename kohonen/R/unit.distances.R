### Calculate distances in a Kohonen map. Crude and
### slow implementation, but hey.

unit.distances <- function(grid, toroidal)
{
  if (!toroidal) {
    if (grid$n.hood == "circular")
      return(as.matrix(stats::dist(grid$pts)))
    if (grid$n.hood == "square")
      return(as.matrix(stats::dist(grid$pts, method="maximum")))
  }

  np <- nrow(grid$pts)
  maxdiffx <- grid$xdim/2
  maxdiffy <- max(grid$pts[,2])/2
  
  result <- matrix(0, np, np)
  for (i in 1:(np-1)) {
    for (j in (i+1):np) {
      diffs <- abs(grid$pts[j,] - grid$pts[i,])
      if (diffs[1] > maxdiffx)
        diffs[1] <- 2*maxdiffx - diffs[1]
      if (diffs[2] > maxdiffy)
        diffs[2] <- 2*maxdiffy - diffs[2]
      
        if (grid$topo == "hexagonal") {
          result[i,j] <- sum(diffs^2)
        } else {
          result[i,j] <- max(diffs)
        }
    }
  }

  if (grid$n.hood == "circular") {
    sqrt(result + t(result))
  } else {
    result + t(result)
  }
}
