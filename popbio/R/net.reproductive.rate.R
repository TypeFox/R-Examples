net.reproductive.rate<-function(A, ...)
{
   if(!is.matrix(A)){stop("A projection matrix is required")}
   A1<-splitA(A, ...)
   Tmat<-A1[[1]]
   Fmat<-A1[[2]]
   s <- length(diag(Tmat))
   N <- try(solve(diag(s) - Tmat), silent=TRUE)
   if(class(N)=="try-error"){r<-NA}
   else{
     R <- Fmat %*% N
     r<-lambda(R)
   }
   r
}
