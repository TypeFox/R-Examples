#' eigen2
#' 
#' Compute eigenvalues and, optionally, eigenvectors of a real symmetric matrix
#' by seraching over ranges of values or ranges of indices.
#' 
#' This new method computes selected eigenvalues and, optionally, eigenvectors
#' of a real symmetric matrix. Eigenvalues and eigenvectors can be selected by
#' specifying either a range of values or a range of indices for the desired
#' eigenvalues.
#' 
#' @param x 
#' symmetric, numeric ddmatrix.
#' @param range 
#' A set of interval endpoints, i.e. a numeric pair.  Controls the
#' set of values over which the eigenvalue search occurs.
#' @param range.type 
#' Controls whether interval \code{range} refers to a set of
#' possible values for the eigenvalues, or a set of indices for the
#' eigenvalues.  Options are "interval" and "index".
#' @param only.values 
#' logical. Determines whether only the eigenvalues should
#' be computed, or if the eigenvectors should as well.
#' @param abstol 
#' The absolute error tolerance for the eigenvalues.
#' @param orfac 
#' Specifies which eigenvectors should be reorthogonalized.
#' Eigenvectors that correspond to eigenvalues which are within
#' tol=orfac*norm(A)of each other are to be reorthogonalized.
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
#' init.grid()
#' 
#' comm.set.seed(seed=1234, diff=TRUE)
#' 
#' x <- crossprod(ddmatrix("rnorm", 10, 3, bldim=2))
#' y <- as.matrix(x)
#' 
#' comm.print(eigen(y))
#' 
#' ### Look for eigenvalues in the range 0 to 10
#' ev <- eigen2(x, range=c(0, 10), only.values=TRUE)
#' comm.print(ev)
#' 
#' finalize()
#' }
#' 
#' @keywords Methods Linear Algebra
#' @name eigen2
#' @export
eigen2 <- function(x, range=c(-Inf, Inf), range.type="interval", only.values=FALSE, abstol=1e-8, orfac=1e-3)
{
    # Basic checking
    must.be(x, "ddmatrix")
    must.be(range, "numeric")
    must.be(range.type, "character")
    must.be(only.values, "logical")
 
    if (x@bldim[1L] != x@bldim[2L])
        comm.stop("The blocking factor for argument 'x' must be square; consider using the redistribute() function")

    # Return eigenvectors or not
    if (only.values)
        jobz <- 'N'
    else
        jobz <- 'V'
    
    # Eigenvalue search
    range.type <- pbdMPI::comm.match.arg(tolower(range.type), c("interval", "index"))
    
    if (range.type == "interval")
    {
        vl <- range[1L]
        vu <- range[2L]
        
        if (vl == -Inf && vu == Inf)
            range <- 'A'
        else
            range <- 'V'
        
        il <- iu <- 0
    }
    else
    {
        il <- range[1L]
        iu <- range[2L]
        
        if (il == -Inf && iu == Inf)
            range <- 'A'
        else
            range <- 'I'
        
        vl <- vu <- 0
    }
    
    
    desca <- base.descinit(dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
    
    out <- base.rpdsyevx(jobz=jobz, range=range, n=x@dim[1L], a=x@Data, desca=desca, vl=vl, vu=vu, il=il, iu=iu, abstol=abstol, orfac=orfac)
    
    if (jobz == 'N')
        return( out )
    else
    {
        if (out$m == 0)
            return( NULL )
        
        z <- new("ddmatrix", Data=out$vectors, dim=x@dim, bldim=x@bldim, ldim=x@ldim, ICTXT=x@ICTXT)
        z <- z[, 1:out$m]
        ret <- list(values=out$values, vectors=z)
        
        return( ret )
    }
}


