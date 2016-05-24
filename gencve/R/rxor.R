rxor <- function(n=1, p=0) {
  y  <- factor(rep(c(1,1,0,0), n))
  x1 <- rep(c(1,0,1,0), n)
  x2 <- rep(c(0,1,1,0), n)
  X <- NULL
  if (p > 0) {
    X <- matrix(ifelse(runif(p*4*n)<0.5, 0, 1), ncol=p, nrow=4*n)
    colnames(X) <- paste0("x", 3:(p+2))
  }
  if (p == 0) {
    out <- data.frame(x1=x1, x2=x2, y=y)
  } else {
    out <- data.frame(x1=x1, x2=x2, X, y=y)
  }
  out
}
