generation.time<-function(A, ...)
{
   if(!is.matrix(A)){stop("A projection matrix is required")}
   A1<-splitA(A, ...)
   Tmat<-A1[[1]]
   Fmat<-A1[[2]]
  
  # probably could add some other checks here
  s <- length(diag(Tmat))
  # check if matrix is singular
   N <- try(solve(diag(s) - Tmat), silent=TRUE)
  if(class(N)=="try-error"){generation.time<-NA}
   else{
     R <- Fmat %*% N
     Ro<- lambda(R)
   lambda <- lambda(A) 
   generation.time = log(Ro)/log(lambda)   
  }
   generation.time
}
