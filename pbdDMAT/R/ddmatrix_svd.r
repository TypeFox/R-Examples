#' Singular Value Decomposition
#' 
#' SVD for distributed matrices with R-like syntax, with
#' calculations performed by the PBLAS and ScaLAPACK libraries.
#' 
#' Extensions of R linear algebra functions.
#' 
#' @param x
#' numeric distributed matrices.
#' @param nu 
#' number of left singular vectors to return when calculating
#' singular values.
#' @param nv 
#' number of right singular vectors to return when calculating
#' singular values.
#' @param LINPACK
#' Ignored.
#' 
#' @return 
#' \code{La.svd()} performs singular value decomposition, and returns the
#' transpose of right singular vectors if any are requested. Singular values
#' are stored as a global R vector. Left and right singular vectors are unique
#' up to sign. Sometimes core R (via LAPACK) and ScaLAPACK will disagree as to
#' what the left/right singular vectors are, but the disagreement is always
#' only up to sign.
#' 
#' \code{svd()} performs singular value decomposition. Differs from
#' \code{La.svd()} in that the right singular vectors, if requested, are
#' returned non-transposed. Singular values are stored as a global R vector.
#' Sometimes core R (via LAPACK) and ScaLAPACK will disagree as to what the
#' left/right singular vectors are, but the disagreement is always only up to
#' sign.
#' 
#' @examples
#' \dontrun{
#' # Save code in a file "demo.r" and run with 2 processors by
#' # > mpiexec -np 2 Rscript demo.r
#' 
#' library(pbdDMAT, quiet = TRUE)
#' init.grid()
#' 
#' # don't do this in production code
#' x <- matrix(1:9, 3)
#' x <- as.ddmatrix(x)
#' 
#' y <- svd(A)
#' print(y)
#' 
#' finalize()
#' }
#' 
#' @aliases svd La.svd
#' @name ddmatrix-svd
#' @rdname ddmatrix-svd
NULL




dmat.svd <- function(x, nu, nv, inplace=FALSE)
{
  ICTXT <- x@ICTXT
  
  # Matrix descriptors
  m <- x@dim[1L]
  n <- x@dim[2L]
  size <- min(m, n)
  bldim <- x@bldim
  
  desca <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=ICTXT)
  
  # U
  if (nu==0){
    jobu <- 'N'
    udim <- c(1L, 1L)
  }
  else {
    jobu <- 'V'
    udim <- c(m, size)
  }
  
  uldim <- base.numroc(dim=udim, bldim=bldim, ICTXT=ICTXT)
  descu <- base.descinit(dim=udim, bldim=bldim, ldim=uldim, ICTXT=ICTXT)
  
  
  # V^T
  if (nv==0){
    jobvt <- 'N'
    vtdim <- c(1L, 1L)
  }
  else {
    jobvt <- 'V'
    vtdim <- c(size, n)
  }
  
  vtldim <- base.numroc(dim=vtdim, bldim=bldim, ICTXT=ICTXT)
  descvt <- base.descinit(dim=vtdim, bldim=bldim, ldim=vtldim, ICTXT=ICTXT)
  
  # Compute 
  out <- base.rpdgesvd(jobu=jobu, jobvt=jobvt, m=m, n=n, a=x@Data, desca=desca, descu=descu, descvt=descvt, inplace=inplace)
  
  if (nu==0)
    u <- NULL
  else {
    u <- new("ddmatrix", Data=out$u, dim=udim, ldim=uldim, bldim=bldim, ICTXT=ICTXT)
    if (nu < u@dim[2L])
      u <- u[, 1L:nu]
  }
  
  if (nv==0)
    vt <- NULL
  else {
    vt <- new("ddmatrix", Data=out$vt, dim=vtdim, ldim=vtldim, bldim=bldim, ICTXT=ICTXT)
    if (nv < vt@dim[1L])
      vt <- vt[1L:nv, ]
  }
  
  ret <- list( d=out$d, u=u, vt=vt )
  
  return( ret )
}







# Games to satisfy codetools' global variable checking --- they don't always play nice with S4
setGeneric(name="La.svd", 
  function(x, ...)
    standardGeneric("La.svd"),
  package="pbdDMAT"
)

setGeneric(name="svd", 
  function(x, ...)
    standardGeneric("svd"),
  package="pbdDMAT"
)




#' @rdname ddmatrix-svd
#' @export
setMethod("La.svd", signature(x="ANY"), 
  function(x, nu = min(n, p), nv = min(n, p)){
    n <- nrow(x)
    p <- ncol(x)
    
    base::La.svd(x=x, nu=nu, nv=nv)
  }
)

#' @rdname ddmatrix-svd
#' @export
setMethod("La.svd", signature(x="ddmatrix"), 
  function(x, nu=min(n, p), nv=min(n, p)) #, ..., inplace=FALSE)
  {
    n <- nrow(x)
    p <- ncol(x)
    
    ret <- dmat.svd(x=x, nu=nu, nv=nv, inplace=FALSE)
    
    return( ret )
  }
)



#' @rdname ddmatrix-svd
#' @export
setMethod("svd", signature(x="ANY"), 
  function(x, nu = min(n, p), nv = min(n, p), LINPACK = FALSE)
  {
    n <- nrow(x)
    p <- ncol(x)
    
    base::svd(x=x, nu=nu, nv=nv, LINPACK=LINPACK)
  }
)

#' @rdname ddmatrix-svd
#' @export
setMethod("svd", signature(x="ddmatrix"), 
  function(x, nu=min(n, p), nv=min(n, p)) #, ..., inplace=FALSE)
  {
    n <- nrow(x)
    p <- ncol(x)
    
    ret <- dmat.svd(x=x, nu=nu, nv=nv, inplace=FALSE)
    
    if (is.ddmatrix(ret$vt))
      ret$vt <- t(ret$vt)
    names(ret)[3] <- "v"
    
    return( ret )
  }
)

