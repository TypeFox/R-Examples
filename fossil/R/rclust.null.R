'rclust.null' <- 
function(groups, dist) {
  c <- length(table(groups))
  dist <- as.matrix(dist)
  #averages table
  atab <- matrix(,10000, c)
  for (i in 1:10000) {
    ng <- sample(groups)
    for (j in 1:c) {
      #for the resampled groups, which equal this 'j'
      z <- which(ng==j)
      atab[i,j] <- sum(dist[z,z])/(length(z)^2-length(z))
    }
  }
  #the final matrix to return, with average win group numbers and sd
  rtmat <- matrix(,c,2)
  for (i in 1:c) rtmat[i,1] <- mean(atab[,i])
  for (i in 1:c) rtmat[i,2] <- sd(atab[,i])
  return(rtmat)
}
