FG <- function(X, weight=NULL, init=NULL, maxiter=100, eps=1e-06, na.action = na.fail)
{
  X <- na.action(X)
  dim.X <- dim(X)
    
  if (length(dim.X)==2) type <- "Matrix"
  if (length(dim.X)==3) type <- "Array"
  if ((length(dim.X) %in% c(2,3))==FALSE) stop("'X' must have two or three dimensions")
    
  if (type == "Array")
     {
      p <- dim.X[1]
      K <- dim.X[3]
      if (dim.X[1] != dim.X[2]) stop("'X' must be an array with dim of the form c(p,p,K)")
      Xt <- aperm(X, c(1,3,2))
      X <- matrix(Xt, ncol=p)
      dim.X <- dim(X)
     }
   
   p <- dim.X[2]
   K <- dim.X[1]/p
   if (floor(K) != ceiling(K)) stop("'X' must be a matrix of k stacked pxp matrices")
  
   if(is.null(weight)) weight <- rep(1,K) 

   if(is.null(init)) init <- diag(p)
   B <- init
    
   res <- .C("FG", as.double(as.vector(X)), as.double(as.vector(B)), as.integer(c(K,p,maxiter)), as.double(as.vector(weight)), as.double(eps), res=double(p^2+1), 
PACKAGE="JADE")$res 

   if(res[p^2+1]==1) stop("try another orthogonal matrix as initial value")
   iter <- res[p^2+1]
   if (iter>maxiter) stop("maxiter reached without convergence")

   V <- matrix(res[1:p^2],p,p)
   D <- X
   for(k in 1:K){
    D[((k-1)*p+1):(k*p),] <- crossprod(V,crossprod(X[((k-1)*p+1):(k*p),],V))
   }
   if(type == "Array"){
    D <- array(t(D),c(p,p,K))
   }
   V
   list(V=V, D=D, iter=iter)
}


