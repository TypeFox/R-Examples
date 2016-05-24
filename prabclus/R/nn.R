"nn" <-
function(distmat,ne=1){
  nc <- ncol(distmat)
  nnd <- c()
  for (i in 1:nc)
    nnd[i] <- sort(distmat[i,])[ne+1]
  out <- mean(nnd)
  out
}
