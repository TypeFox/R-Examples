fundamental.matrix<-function(A, ...)
{
   if(!is.matrix(A)){stop("A projection matrix is required")}
   A1<-splitA(A, ...)
   Tmat<-A1[[1]]
 
   s <- length(diag(Tmat))
   # check if matrix is singular 
   N <- try(solve(diag(s) - Tmat), silent=TRUE)
  if(class(N)=="try-error"){fundamental.matrix<-"Transition matrix is singular"}
   else{

   
   var<- (2*diag(diag(N)) - diag(s)) %*% N - N * N
     dimnames(var)<-dimnames(A)
   total <- margin.table(N,2)
   vareta <- margin.table(2*(N %*% N) - N, 2) - total * total

   fundamental.matrix <- list(
     N = N,
     var = var,
     cv = var^.5/N,
     meaneta = total,
     vareta = vareta
                )
   }
   fundamental.matrix
}
