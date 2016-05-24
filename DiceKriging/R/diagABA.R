diagABA <- function(A, B){
  D <- B
  D[lower.tri(D, diag=FALSE)] <- 0
  D[upper.tri(D, diag=FALSE)] <- 2*D[upper.tri(D, diag=FALSE)]
  D <- (A%*%D)*A
  rowSums(D)
}