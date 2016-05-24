
uniclogit <- function(g, Y, w, n, group){
  theta <- rep(0.1, 1000)
  k <- 2
  firstdev <- 1
  while (abs(firstdev) > 1e-5) {
    firstdev <- score(g, Y, w, n, group, theta[k - 1])
    theta[k] <- theta[k - 1] - firstdev / seconddev(g, Y, w, n, group, theta[k - 1])
    k <- k + 1
  }
  return(theta[max(which(theta != 0.1))])
}