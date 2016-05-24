sensitivity<-function(A, zero=FALSE)
{
   ev <- eigen(A)
   lmax <- which.max(Re(ev$values))
   W <- ev$vectors
   w <- abs(Re(W[, lmax]))
   V <- try(Conj(solve(W)), silent=TRUE)
   if (class(V) == "try-error"){stop("matrix A is singular") }
   v <- abs(Re(V[lmax, ]))
   s <- v %o% w
   if (zero) {
        s[A == 0] <- 0
   }
   dimnames(s) <- dimnames(A)
   s
}
