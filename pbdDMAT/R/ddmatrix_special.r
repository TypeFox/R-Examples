#' Generate Companion Matrices
#' 
#' Methods for constructing companion matrices of an n-degree polynomial.
#' 
#' For a degree n polynomial,
#' 
#' \eqn{x^n + a_{n-1}x^{n-1} + \dots + a_1x + a_0}
#' 
#' its associated companion matrix is a matrix of the form
#' 
#' \eqn{\left[\begin{array}{cccccc} 0 & 0 & 0 & \dots & 0 & -a_0\\ 1 & 0 & 0 &
#' \dots & 0 & -a_1\\ 0 & 1 & 0 & \dots & 0 & -a_2\\ \vdots & \vdots & \vdots &
#' \ddots & \vdots & \vdots\\ 0 & 0 & 0 & \dots & 1 & -a_{n-1}
#' \end{array}\right]}
#' 
#' In the function call, we assume that the argument '\code{coef}' is ordered
#' from \eqn{a_0} to \eqn{a_{n-1}}.
#' 
#' NOTE that we assume that the leading coefficient is 1.
#' 
#' @param coef 
#' Vector of polynomial coefficients, listed in increasing order
#' (by index; see details below).
#' @param type 
#' "matrix" or "ddmatrix".
#' @param ... 
#' Additional arguments.
#' @param bldim 
#' blocking dimension.
#' @param ICTXT 
#' BLACS context number.
#' 
#' @return 
#' Returns a matrix or a distributed matrix.
#' 
#' @keywords Data Generation
#' 
#' @export
companion <- function(coef, type="matrix", ..., bldim=.pbd_env$BLDIM, ICTXT=.pbd_env$ICTXT)
{
  type <- pbdMPI::comm.match.arg(type, c("matrix", "ddmatrix"))
  
  if (type=="ddmatrix"){
    if (length(bldim)==1)
      bldim <- rep(bldim, 2L)
    
    dim <- rep(length(coef), 2L)
    ldim <- base.numroc(dim=dim, bldim=bldim, ICTXT=ICTXT)
    
    descx <- base.descinit(dim=dim, bldim=bldim, ldim=ldim, ICTXT=ICTXT)
    
    out <- base.pdmkcpn1(coef=coef, descx=descx)
    ret <- new("ddmatrix", Data=out, dim=dim, ldim=ldim, bldim=bldim, ICTXT=ICTXT)
  }
  else
  {
    n <- length(coef)
    ret <- cbind(rbind(rep(0, n-1), diag(1, nrow=n-1, ncol=n-1)), -coef)
  }
  
  return( ret )
}



#' Generate Hilbert Matrices
#' 
#' Methods for constructing Hilbert matrices: H[i,j] = 1/(i+j-1)
#' 
#' This constructs the square Hilbert matrix of order \code{n}. The return is
#' either a matrix or a distributed matrix depending on the argument
#' \code{type=}.
#' 
#' @param n 
#' number of rows and columns.
#' @param type 
#' "matrix" or "ddmatrix".
#' @param ... 
#' Additional arguments.
#' @param bldim 
#' blocking dimension.
#' @param ICTXT 
#' BLACS context number.
#' 
#' @return 
#' Returns a matrix or a distributed matrix.
#' 
#' @keywords Data Generation
#' 
#' @examples
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' dx <- Hilbert(100, type="ddmatrix")
#' 
#' print(dx)
#' 
#' finalize()
#' }
#' 
#' @export
Hilbert <- function(n, type="matrix", ..., bldim=.pbd_env$BLDIM, ICTXT=.pbd_env$ICTXT)
{
  type <- pbdMPI::comm.match.arg(type, c("matrix", "ddmatrix"))
  if (type == "ddmatrix")
  {
    dim <- c(n, n)
    ldim <- base.numroc(dim=dim, bldim=bldim, ICTXT=ICTXT, fixme=TRUE)
    descx <- base.descinit(dim=dim, bldim=bldim, ldim=ldim, ICTXT=ICTXT)
    
    out <- base.pdhilbmk(descx)
    
    ret <- new("ddmatrix", Data=out, dim=dim, ldim=ldim, bldim=bldim, ICTXT=ICTXT)
  }
  else
  {
    ret <- base.dhilbmk(n)
  }
  
  return( ret )
}
