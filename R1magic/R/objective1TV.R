objective1TV <- function(x,T,phi,y,lambda) {
  # Part of R1Magic by mehmet.suzen@physics.org
  X  <- T %*% x
  OF <- norm( phi %*% X - y , c("F") ) ^ 2 +  lambda * TV1 (X)
 return (OF)
}
