`R11.AR1` <-
function(lambda0,dimension) {
  temp <- matrix(NA,dimension,dimension)
  diag(temp) <- 1
  for (i in 2:dimension) {
    for (j in 1:(i-1)) {
      temp[i,j] <- lambda0^(i-j)
      temp[j,i] <- temp[i,j]
    }
  }
  return(temp)
}

