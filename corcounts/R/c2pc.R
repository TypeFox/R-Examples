`c2pc` <-
function(Cin) {
  ntemp <- dim(Cin)[2]
  upper.triangle.index <- matrix(0,(ntemp-1)*ntemp/2,ntemp)
  counter <- 1
  for (i in 1:(ntemp-1)) {
    for (j in (i+1):ntemp) {
      upper.triangle.index[counter,1:2] <- c(i,j)
      if (i>1) { upper.triangle.index[counter, 3:(1+i)] <- 1:(i-1) }
      counter <- counter + 1
    }
  }

  upper.triangle.values <- double((ntemp-1)*ntemp/2)
  for (k in 1:(ntemp-1)) {
    upper.triangle.values[k] <- Cin[upper.triangle.index[k,1],upper.triangle.index[k,2]]
  }

  if (ntemp<(dim(upper.triangle.index)[1]+1)) {
    for (k in ntemp:dim(upper.triangle.index)[1]) {
      upper.triangle.values[k] <- berechne.partial.corr(k,Cin,ntemp,upper.triangle.index,upper.triangle.values)
    }
  }

  Theta <- matrix(NA,ntemp,ntemp)
  counter <- 1
  for (i in 2:ntemp) {
    for (j in (i-1):(ntemp-1)) {
      Theta[i-1,j+1] <- upper.triangle.values[counter]
      counter <- counter + 1
    }
  }

  return(Theta)
}

