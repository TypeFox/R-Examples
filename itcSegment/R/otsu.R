
otsu <- function(y,m=1){

  yvals <- sort(unique(y))
  L <- length(yvals)
  per <- as.vector(table(y)) / length(y)

  P <- matrix(0, nrow=L, ncol=L)
  S <- matrix(0, nrow=L, ncol=L)
  H <- matrix(0, nrow=L, ncol=L)

  P[1,] <- cumsum(per)
  S[1,] <- cumsum(per * yvals[1:L])
  for(u in 2:L)
    for(v in u:L){
      P[u,v] <- P[1,v] - P[1,u-1]
      S[u,v] <- S[1,v] - S[1,u-1]
    }
  H <- S^2 / P


  x <- seq(L)
  n <- length(x)
  if (n -1 < m)
    stop("The number of thresholds is larger than the unique values minus 1.")
  e <- 0
  h <- m
  a <- 1:m
  rule <- c(0, a, L)
  sigma2 <- sum(sapply(1:(m+1), function(i) H[rule[i]+1, rule[i+1]]))

  thresh <- yvals[a]
  nmmp1 <- n - m + 1
  mp1 <- m + 1
  while (a[1] != nmmp1) {
    if (e < n - h) {
      h <- 1
      e <- a[m]
      j <- 1
    }
    else {
      h <- h + 1
      e <- a[mp1 - h]
      j <- 1:h
    }
    a[m - h + j] <- e + j
    if(a[m] != L){
      rule <- c(0, a, L)
      new <- sum(sapply(1:(m+1), function(i) H[rule[i]+1, rule[i+1]]))
      if(new > sigma2){
        sigma2 <- new
        thresh <- yvals[a]
      }
    }
  }
  thresh
}
