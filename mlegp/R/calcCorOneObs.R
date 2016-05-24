`calcCorOneObs` <-
function(Z, B, a, newZ) {
  obs = dim(Z)[1]
  newCor = matrix(0, ncol = obs)
  for (i in 1:obs) {
        z <- as.numeric(Z[i,]) 
        newCor[i] <- exp(-sum(B*abs(z-newZ)**a))
  }
  return (newCor)
}

