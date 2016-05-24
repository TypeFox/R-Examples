# ------------------------------------------------
# Linear Equations
# ------------------------------------------------

#' rpdgetri
#' 
#' Matrix inversion.
#' 
#' For advanced users only.
#' 
#' @param n
#' Problem size.
#' @param a
#'  Matrix.
#' @param desca
#' ScaLAPACK descriptor array.
#' 
#' @export
base.rpdgetri <- function(n, a, desca)
{
    aldim <- base.numroc(desca[3:4], desca[5:6], ICTXT=desca[2])
    
    if (!is.double(a))
        storage.mode(a) <- "double"
    
    out <- .Call(R_PDGETRI, a, as.integer(desca))
    
    if (out$info!=0)
        pbdMPI::comm.warning(paste("ScaLAPACK returned INFO=", out$info, "; returned solution is likely invalid", sep=""))
    
    return( out$A )
}



#' rpdgesv
#' 
#' Solving a (square) system of equations.
#' 
#' For advanced users only.
#' 
#' @param n
#' Problem size.
#' @param nrhs
#' Number of right hand sides.
#' @param a,b
#' Matrix.
#' @param desca,descb
#' ScaLAPACK descriptor array.
#' 
#' @export
base.rpdgesv <- function(n, nrhs, a, desca, b, descb)
{
    aldim <- base.numroc(desca[3:4], desca[5:6], ICTXT=desca[2])
    bldim <- base.numroc(descb[3:4], descb[5:6], ICTXT=descb[2])
    
    # max of the local dimensions
    mxldims <- c(base.maxdim(aldim), base.maxdim(bldim))
    
    if (!is.double(a))
        storage.mode(a) <- "double"
    if (!is.double(b))
        storage.mode(b) <- "double"
    
    # Call ScaLAPACK
    out <- .Call(R_PDGESV,
                 as.integer(n), as.integer(nrhs), as.integer(mxldims),
                 a, as.integer(desca), b, as.integer(descb))
    
    if (out$info!=0)
        pbdMPI::comm.warning(paste("ScaLAPACK returned INFO=", out$info, "; returned solution is likely invalid", sep=""))
    
    return( out$B ) 
}



# ------------------------------------------------
# Matrix Factorizations
# ------------------------------------------------

#' rpdgesvd
#' 
#' SVD.
#' 
#' For advanced users only.
#' 
#' @param jobu,jobvt
#' Control for u/vt return.
#' @param m,n
#' Problem size.
#' @param a
#'  Matrix.
#' @param desca,descu,descvt
#' ScaLAPACK descriptor array.
#' @param ...
#' Ignored
#' @param inplace
#' Should the computation be done in-place or not.  For REALLY advanced users only.
#' 
#' @export
base.rpdgesvd <- function(jobu, jobvt, m, n, a, desca, descu, descvt, ..., inplace=FALSE)
{
    size <- min(m, n)
    
    aldim <- dim(a)
    uldim <- base.numroc(descu[3:4], descu[5:6], ICTXT=descu[2])
    vtldim <- base.numroc(descvt[3:4], descvt[5:6], ICTXT=descvt[2])
    
    mxa <- pbdMPI::allreduce(max(aldim), op='max')
    mxu <- pbdMPI::allreduce(max(uldim), op='max')
    mxvt <- pbdMPI::allreduce(max(vtldim), op='max')
    
    if (all(aldim==1))
        desca[9L] <- mxa
    if (all(uldim==1))
        descu[9L] <- mxu
    if (all(vtldim==1))
        descvt[9L] <- mxvt
    
    if (desca[3L]>1){
        if (pbdMPI::allreduce(desca[9L], op='max')==1)
            desca[9L] <- mxa
    }
    if (descu[3L]>1){
        if (pbdMPI::allreduce(descu[9L], op='max')==1)
            desca[9L] <- mxu
    }
    if (descvt[3L]>1){
        if (pbdMPI::allreduce(descvt[9L], op='max')==1)
            desca[9L] <- mxvt
    }
    
    if (!is.double(a))
        storage.mode(a) <- "double"
    
    # FIXME currently does nothing
    if (inplace)
        inplace <- 'Y'
    else
        inplace <- 'N'
    
    # Call ScaLAPACK
    out <- .Call(R_PDGESVD, 
                        as.integer(m), as.integer(n), as.integer(size),
                        a, as.integer(desca), 
                        as.integer(uldim), as.integer(descu),
                        as.integer(vtldim), as.integer(descvt),
                        as.character(jobu), as.character(jobvt), inplace)
    
    if (out$info!=0)
        pbdMPI::comm.warning(paste("ScaLAPACK returned INFO=", out$info, "; returned solution is likely invalid", sep=""))
    
    ret <- list( d=out$d, u=out$u, vt=out$vt )
    
    return( ret )
}



#' rpdsyev
#' 
#' Symmetric eigenvalue decomposition.
#' 
#' For advanced users only.
#' 
#' @param jobz
#' Control for if vectors/values/both are returned.
#' @param uplo
#' Triangle where the information is stored (in the symmetric matrix).
#' @param n
#' Problem size.
#' @param a
#' Matrix.
#' @param desca,descz
#' ScaLAPACK descriptor array.
#' 
#' @export
base.rpdsyev <- function(jobz, uplo, n, a, desca, descz)
{
    aldim <- dim(a)
    zldim <- base.numroc(descz[3:4], descz[5:6], ICTXT=descz[2])
    
    mxa <- pbdMPI::allreduce(max(aldim), op='max')
    mxz <- pbdMPI::allreduce(max(zldim), op='max')
    
    if (all(aldim==1))
        desca[9L] <- mxa
    if (all(zldim==1))
        descz[9L] <- mxz
    
    if (!is.double(a))
        storage.mode(a) <- "double"
    
    # Call ScaLAPACK
    out <- .Call(R_PDSYEV, 
                        as.character(jobz), as.character(uplo),
                        as.integer(n), a, as.integer(desca), as.integer(aldim),
                        as.integer(zldim), as.integer(descz))
    
    if (out$info!=0)
        pbdMPI::comm.warning(paste("ScaLAPACK returned INFO=", out$info, "; returned solution is likely invalid", sep=""))
    
    out$values <- rev(out$values)
    out$info <- NULL
    
    return( out )
}



#' rpdpotrf
#' 
#' Cholesky factorization.
#' 
#' For advanced users only.
#' 
#' @param uplo
#' Triangle where the information is stored (in the symmetric matrix).
#' @param n
#' Problem size.
#' @param a
#'  Matrix.
#' @param desca
#' ScaLAPACK descriptor array.
#' 
#' @export
base.rpdpotrf <- function(uplo, n, a, desca)
{
    
    if (!is.double(a))
        storage.mode(a) <- "double"
    
    # Call ScaLAPACK
    ret <- .Call(R_PDPOTRF,
                 as.integer(n), a, as.integer(desca),
                 as.character(uplo))
    
    if (ret$info!=0)
        pbdMPI::comm.warning(paste("ScaLAPACK returned INFO=", ret$info, "; returned solution is likely invalid", sep=""))
    
    return( ret ) 
}



#' rpdsyevx
#' 
#' Genearlized eigenvalue problem.
#' 
#' For advanced users only.
#' 
#' @param jobz
#' Control for if vectors/values/both are returned.
#' @param range
#' Parameter to determine the search criteria for eigenvalues.
#' @param n
#' Problem size.
#' @param a
#'  Matrix.
#' @param desca
#' ScaLAPACK descriptor array.
#' @param vl,vu
#' Endpoints of the interval subset of the real line in which to search for eigenvalues, if specified by \code{range}.
#' @param il,iu
#' Eigenvalues with indices \code{il}, ..., \code{iu} will be found, if specified by \code{range}.
#' @param abstol
#' Absolute error tolerance for the eigenvalues.
#' @param orfac
#' Eigenvectors with eigenvalues below orfac*norm(a) of each other are reorthogonalized.
#' 
#' @export
base.rpdsyevx <- function(jobz, range, n, a, desca, vl, vu, il, iu, abstol=1e-8, orfac=1e-3)
{
    
    if (!is.double(a))
        storage.mode(a) <- "double"
    
    
    # Call ScaLAPACK
    ret <- .Call(R_PDSYEVX,
               as.character(jobz), as.character(range),
               as.integer(n), a, as.integer(desca),
               as.double(vl), as.double(vu), as.integer(il), as.integer(iu),
               as.double(abstol), as.double(orfac))
    
    return( ret )
}


#' rpdgetrf
#' 
#' LU factorization.
#' 
#' For advanced users only.
#' 
#' @param a
#'  Matrix.
#' @param desca
#' ScaLAPACK descriptor array.
#' 
#' @export
base.rpdgetrf <- function(a, desca)
{
    m <- desca[3L]
    n <- desca[4L]
    
    aldim <- dim(a)
    
    mxrow <- base.igamx2d(ICTXT=desca[2L], SCOPE='Row', m=1L, n=1L, x=aldim[1L], lda=1L, RDEST=-1L, CDEST=-1L)
    lipiv <- mxrow + desca[5L]
#    lipiv <- base.maxdim(aldim)[1L] + desca[5L]
    
    if (!is.double(a))
        storage.mode(a) <- "double"
    
    # Call ScaLAPACK
    out <- .Call(R_PDGETRF,
                 as.integer(m), as.integer(n),
                 a, as.integer(aldim), as.integer(desca),
                 as.integer(lipiv))
    
    if (out$info!=0)
        pbdMPI::comm.warning(paste("ScaLAPACK returned INFO=", out$info, "; returned solution is likely invalid", sep=""))
    
    return( out$A ) 
}



# ------------------------------------------------
# Auxillary
# ------------------------------------------------

#' rpdlange
#' 
#' Matrix norms.
#' 
#' For advanced users only.
#' 
#' @param norm
#' Type of norm.
#' @param m,n
#' Problem size
#' @param a
#' Matrix.
#' @param desca
#' ScaLAPACK descriptor array.
#' 
#' @export
base.rpdlange <- function(norm, m, n, a, desca)
{
    if (length(norm)>1L)
        norm <- norm[1L]
    
    norm <- toupper(norm)
    
    if (!is.double(a))
        storage.mode(a) <- "double"
    
    ret <- .Call(R_PDLANGE, 
                norm, as.integer(m), as.integer(n),
                a, as.integer(desca))
    
    return( ret )
}



#' rpdtrcon
#' 
#' Inverse condition number of a triangular matrix.
#' 
#' For advanced users only.
#' 
#' @param norm
#' Type of norm.
#' @param uplo
#' Triangle where information is stored.
#' @param diag
#' Specifies if the matrix is unit triangular or not.
#' @param n
#' Problem size
#' @param a
#' Matrix.
#' @param desca
#' ScaLAPACK descriptor array.
#' 
#' @export
base.rpdtrcon <- function(norm, uplo, diag, n, a, desca)
{
    if (length(norm)>1L)
        norm <- norm[1L]
    
    norm <- toupper(norm)
    uplo <- toupper(uplo)
    diag <- toupper(diag)
    
    if (!is.double(a))
        storage.mode(a) <- "double"
    
    ret <- .Call(R_PDTRCON, 
                norm, uplo, diag, 
                as.integer(n), a, as.integer(desca))
    
    if (ret[2L] < 0)
        pbdMPI::comm.warning(paste("INFO =", ret[2L]))
    
    return( ret[1L] )
}



#' rpdgecon
#' 
#' Inverse condition number of a general matrix.
#' 
#' For advanced users only.
#' 
#' @param norm
#' Type of norm.
#' @param m,n
#' Problem size
#' @param a
#' Matrix.
#' @param desca
#' ScaLAPACK descriptor array.
#' 
#' @export
base.rpdgecon <- function(norm, m, n, a, desca)
{
    if (length(norm)>1L)
        norm <- norm[1L]
    
    norm <- toupper(norm)
    
    if (!is.double(a))
        storage.mode(a) <- "double"
    
    ret <- .Call(R_PDGECON, norm, as.integer(m), as.integer(n), a, as.integer(desca))
    
    if (ret[2] < 0)
        pbdMPI::comm.warning(paste("INFO =", ret[2]))
    
    return( ret[1] )
}



# ------------------------------------------------
# Utility
# ------------------------------------------------

#' rpdgemr2d
#' 
#' General 2d block cyclic redistribution function.
#' 
#' For advanced users only.
#' 
#' @param x
#' Matrix.
#' @param descx,descy
#' ScaLAPACK descriptor array.
#' 
#' @export
base.rpdgemr2d <- function(x, descx, descy)
{
    ldimy <- base.numroc(dim=descy[3L:4L], bldim=descy[5L:6L], ICTXT=descy[2L])
    ldimy <- as.integer(ldimy)
    descx <- as.integer(descx)
    descy <- as.integer(descy)
    m <- descx[3L]
    n <- descx[4L]
    
    # context 0 is always passed since pxgemr2d 
    # requires the grids to have at least 1 processor in common
    ### TODO integrate PIGEMR2D
    if (!is.double(x))
        storage.mode(x) <- "double"
    
    if (!base.ownany(dim=c(m, n), bldim=descy[5L:6L], ICTXT=descy[2L]))
        ret <- matrix(0.0, 1L, 1L)
    
    return( ret )
}
