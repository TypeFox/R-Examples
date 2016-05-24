createDesMat <- function(n, p) {
  X0 <- rep(1, n*p) # first row 
  X11 <- kronecker(diag(1,n-1), t(rep(1,p)))
  X12 <- matrix(rep(-1, p*(n-1)), nrow=n-1, ncol=p)
  
  X1 <- cbind(X11, X12) # first block
  
  X21 <- cbind(diag(1, p-1), rep(-1, p-1))
  X22 <- X21
  
  for(i in 1:(n-1)) X21 <- cbind(X21, X22) # second block
  
  rbind(X0, X1, X21)
}