rclust.weights <-
function(groups, dist) {
  c <- length(table(groups))
  D <- as.matrix(dist)
  n <- dim(D)[1]
  weights <- matrix(,n,c)
  for (k in 1:n) {
    for (i in 1:c) {
      z <- which(groups==i)
      if (groups[k]==i) weights[k,i] <- sum(D[k,z])/(length(z)-1)
      else weights[k,i] <- sum(D[k,z])/length(z)
    }
  }
  return(weights)
}
