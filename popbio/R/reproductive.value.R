reproductive.value<-function(A)
{
   ev <- eigen(A)
   lmax <- which.max(Re(ev$values))
   W <- ev$vectors
   w <- abs(Re(W[, lmax]))
   V <- try(Conj(solve(W)), silent=TRUE)
   if (class(V) == "try-error"){stop("matrix A is singular") }
   v <- abs(Re(V[lmax, ]))
   names(v) <-colnames(A)
   v/v[1]
}
