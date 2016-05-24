eiv <- function(L,  identity=seq(NCOL(L)), ...){
   A1 <- L[ identity, , drop=FALSE]
   g <- solve(A1)
   if(1e-14 < max(abs(diag(1, length(identity)) - A1 %*% g)))
      warning("Inverse is not well conditioned. Consider setting identity to select different rows.")
   B <- array(NA, dim(L))
   B[identity, ] <- diag(1, length(identity))
   B[-identity,] <- L[-identity,, drop=FALSE] %*% g
   dimnames(B) <- list(dimnames(L)[[1]], paste("factor", seq(NCOL(L))))
   list(loadings=B, Th=t(A1), method="eiv", orthogonal=FALSE, convergence=TRUE,
        Phi= tcrossprod(A1))
   }
