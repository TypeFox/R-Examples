echelon <- function(L, reference=seq(NCOL(L)), ...) {

   # Split L in reference part and the rest
   A1 <- L[reference,, drop=FALSE]
   #A2 is L[-reference,]

   # Compute the part of A Phi A' corresponding to the reference variables
   # Compute cholesky rot = rotated reference part
   # No check or error message for singularity. Exact singularity is rare in
   # practice but ill-conditioning is a real danger.

   #  now assuming orthogonal (Phi=I)
   #newPhi <- if (is.null(Phi)) A1 %*% t(A1) else A1 %*% Phi %*% t(A1)
   #B1 <- t(chol(newPhi))
   
   B1 <- t(chol(A1 %*% t(A1)))

   # Transformation matrix: B1 = A1 * Tmat
   # Rotated solution for non-reference part: B2 = A2 * Tmat
   Tmat <- solve(A1, B1)
   
   # Assemble rotated solution
   B <- matrix(0, NROW(L), NCOL(L))
   B[reference,]  <- B1
   B[-reference,] <- L[-reference,, drop=FALSE] %*% Tmat

   dimnames(B) <- list(dimnames(L)[[1]], paste("factor", seq(NCOL(L))))
   list(loadings=B, Th=Tmat, method="echelon", orthogonal=TRUE, 
       convergence=TRUE)
}
