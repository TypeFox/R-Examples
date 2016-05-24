#' Matrix Exponentiation
#' 
#' Routines for matrix exponentiation.
#' 
#' Formally, the exponential of a square matrix \code{X} is a power series:
#' 
#' \eqn{expm(X) = id + X/1! + X^2/2! + X^3/3! + \dots}
#' 
#' where the powers on the matrix correspond to matrix-matrix multiplications.
#' 
#' \code{expm()} directly computes the matrix exponential of a distributed,
#' dense matrix.  The implementation uses Pade' approximations and a
#' scaling-and-squaring technique (see references).
#' 
#' @references
#' "New Scaling and Squaring Algorithm for the Matrix Exponential"
#' Awad H. Al-Mohy and Nicholas J. Higham, August 2009
#' 
#' @param x 
#' A numeric matrix or a numeric distributed matrix.
#' @param t
#' Scaling parameter for x.
#' @param p
#' Order of the Pade' approximation.
#' 
#' @return 
#' Returns a distributed matrix.
#' 
#' @examples
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' x <- matrix("rnorm", 5, 5, bldim=2)
#' expm(x)
#' 
#' }
#' @keywords Methods Linear Algebra
#' @name expm
#' @rdname expm
setGeneric(name="expm", 
  function(x, t=1, p=6) 
    standardGeneric("expm"), 
  package="pbdDMAT"
)





p_matpow_by_squaring <- function(A, b=1)
{
  b <- as.integer(b)
  
  desca <- base.descinit(dim=A@dim, bldim=A@bldim, ldim=A@ldim, ICTXT=A@ICTXT)
  
  out <- base.p_matpow_by_squaring_wrap(A=A@Data, desca=desca, b=b)
  
  ret <- new("ddmatrix", Data=out, dim=A@dim, ldim=A@ldim, bldim=A@bldim, ICTXT=A@ICTXT)
  
  return( ret )
}


#matpow <- function(A, n)
#{
#  m <- nrow(A)
#  
##  if (n==0)
##    return(diag(1, m))
##  else if (n==1)
##    return(A)
##  
##  if (n >= 2^8)
##  {
##    E <- eigen(A)
##    B <- E$vectors %*% diag(E$values^n) %*% solve(E$vectors)
##    
##    return( B )
##  }
#  
#  B <- matpow_by_squaring(A, n)
#  
#  return(B)
#}



p_matexp_pade <- function(A, p)
{
  desca <- base.descinit(dim=A@dim, bldim=A@bldim, ldim=A@ldim, ICTXT=A@ICTXT)
  
  out <- base.p_matexp_pade_wrap(A=A@Data, desca=desca, p=p)
  
  N <- new("ddmatrix", Data=out$N, dim=A@dim, ldim=A@ldim, bldim=A@bldim, ICTXT=A@ICTXT)
  D <- new("ddmatrix", Data=out$D, dim=A@dim, ldim=A@ldim, bldim=A@bldim, ICTXT=A@ICTXT)
  
  R <- solve(D) %*% N
  
  return( R )
}



#' @rdname expm
#' @export
setMethod("expm", signature(x="matrix"), 
  function(x, t=1, p=6)
  {
    if (nrow(x) != ncol(x))
      stop("Matrix exponentiation is only defined for square matrices.")
    
    if (p > 13)
      stop("argument 'p' must be between 1 and 13.")
    
    S <- base.matexp(A=x, p=p, t=t)
    
    return( S )
  }
)


## TODO export this from base src/
matexp_scale_factor <- function(x)
{
#  theta <- c(3.7e-8, 5.3e-4, 1.5e-2, 8.5e-2, 2.5e-1, 5.4e-1, 9.5e-1, 1.5e0, 2.1e0, 2.8e0, 3.6e0, 4.5e0, 5.4e0, 6.3e0, 7.3e0, 8.4e0, 9.4e0, 1.1e1, 1.2e1, 1.3e1)
  theta <- c(1.5e-2, 2.5e-1, 9.5e-1, 2.1e0, 5.4e0)
  
  
  # 1-norm
  x_1 <- norm(x, type="O") # max(colSums(abs(x))) 
  
  for (th in theta)
  {
    if (x_1 <= th)
      return( 0 )
  }
  
  j <- ceiling(log2(x_1/theta[5]))
  n <- 2^j
  
  return( n )
}




#' @rdname expm
#' @export
setMethod("expm", signature(x="ddmatrix"), 
  function(x, t=1, p=6)
  {
    if (nrow(x) != ncol(x))
      stop("Matrix exponentiation is only defined for square matrices.")
    
    n <- matexp_scale_factor(x)
    
    if (n == 0)
      return( p_matexp_pade(t*x, p=p) )
    
    x <- t*x/n
    
    S <- p_matexp_pade(x, p=p)
    S <- p_matpow_by_squaring(S, n)
    
    return( S )
  }
)


