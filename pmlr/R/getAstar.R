# computes the Fisher information A* for the penalized likelihood

getAstar <- function(p, J, A, W, V) {

   if (det(A) < 1e-05) diag(A) <- diag(A) + 1e-05

   Ainv <- solve(A)
   T <- W %*% cbind(as.vector(Ainv))
   R1 <- V %*% (cbind(as.vector(Ainv)) %x% diag(p*J)); R2 <- W %*% (Ainv %*% T %x% diag(p*J))

   Astar <- A - 0.5 * (R1 - R2)
   return(Astar)

}
