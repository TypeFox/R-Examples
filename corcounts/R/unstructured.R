`unstructured` <-
function(dimension) {
  temp <- matrix(NA,dimension,dimension)
  diag(temp) <- 1
  for (i in 2:dimension) {
    for (j in 1:(i-1)) {
      temp[i,j] <- runif(1,-.9,.9)
      temp[j,i] <- temp[i,j]
    }
  }
  return(pc2c(temp))
}

