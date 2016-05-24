
geigen <- function(A, B, symmetric, only.values=FALSE) {

    if(!is.matrix(A)) stop("Argument A should be a matrix")
    if(!is.matrix(B)) stop("Argument B should be a matrix")
    dimA <- dim(A)
    if(dimA[1]!=dimA[2]) stop("A must be a square matrix")
    dimB <- dim(B)
    if(dimB[1]!=dimB[2]) stop("B must be a square matrix")
    if(dimA[1]!=dimB[1]) stop("A and B must have the same dimensions")

    if(dimA[1]==0) stop("Matrix A has zero rows/columns")
    if(dimB[1]==0) stop("Matrix B has zero rows/columns")

    if(!all(is.finite(A))) stop("Matrix A may not contain infinite/NaN/NA")
    if(!all(is.finite(B))) stop("Matrix B may not contain infinite/NaN/NA")

    if( missing(symmetric) ) {
        # test if both matrices are symmetric
        # both should be symmetric  for symmetric algorithm
        # expensive test
        symmetric <- isSymmetric(A) && isSymmetric(B)
    }

    complex.A <- is.complex(A)
    complex.B <- is.complex(B)
    complex.AB <- FALSE
    if( complex.A || complex.B ) {
        # at least one matrix is complex ==> use complex routines
        if( !complex.A ) storage.mode(A) <- "complex"
        if( !complex.B ) storage.mode(B) <- "complex"
        complex.AB <- TRUE
    } else {
        # none of the matrices is complex ==> use double routines
        if( !is.double(A) ) storage.mode(A) <- "double"
        if( !is.double(B) ) storage.mode(B) <- "double"
    }

    if( symmetric ) {
        if( complex.AB ) {
            geigen.zhegv(A, B, only.values)
        } else {
            geigen.dsygv(A, B, only.values)
        }
    } else {
        if( complex.AB ) {
            geigen.zggev(A, B, only.values)
        } else {
            geigen.dggev(A, B, only.values)
        }
    }
}

geigen.dggev <- function(A,B, only.values=FALSE) {

    # interface to Lapack dggev
    # for generalized eigenvalue problem
    # general matrices
    # for symmetric matrices use dsygv

    # jobvl.char <- "N"
    # jobvr.char <- if(!only.values) "V" else "N"

    kjobvl <- 1L
    kjobvr <- if(!only.values) 2L else 1L

    dimA <- dim(A)
    n <- dimA[1]

    # calculate optimal workspace
    lwork <- -1L
    work  <- numeric(1)
    z <- .Fortran("xdggev", kjobvl, kjobvr, n, A, B, numeric(1), numeric(1),
                            numeric(1), numeric(1), 1L, numeric(1), n,
                            work=work, lwork,info=integer(1L))
    lwork <- as.integer(z$work[1])
    lwork <- max(lwork, as.integer(8*n))
    work  <- numeric(lwork)

    if( !only.values )
        z <- .Fortran("xdggev", kjobvl, kjobvr, n, A, B, alphar=numeric(n), alphai=numeric(n),
                                beta=numeric(n), numeric(1),1L, vr=matrix(0,nrow=n,ncol=n), n,
                                work, lwork,info=integer(1L), PACKAGE="geigen")
    else
        z <- .Fortran("xdggev", kjobvl, kjobvr, n, A, B, alphar=numeric(n), alphai=numeric(n),
                                beta=numeric(n), numeric(1), 1L,  numeric(1), 1L,
                                work, lwork,info=integer(1L), PACKAGE="geigen")

    if( z$info != 0 ) .ggev_Lapackerror(z$info,n)

    # warning: z$beta may contain exact zeros
    # so dividing by z$beta can result in Inf values

    if( all(z$alphai==0) ) {
        values <- z$alphar/z$beta
        alpha  <- z$alphar
        if(!only.values) vectors <- z$vr
    }
    else {
        # at least one alphai is <> 0
        # alphai == 0 implies real eigenvalue
        # alphai  > 0 implies a complex conjugate pair of eigenvalues
        #             alphai[j] > 0 ==> alphai[j+1] < 0
        
        # use specialized complex division function in order to avoid NaN+NaNi 
        # on Mac OS X with the official R (because of broken gcc 4.2.X)

        alpha  <- complex(real=z$alphar, imaginary=z$alphai) 
        values <- complexdiv(alpha,z$beta)

        if( !only.values ) {
            vectors <- z$vr
            storage.mode(vectors) <- "complex"
            for( j in which(z$alphai>0) ) {
                vectors[,j]   <- complex(real=z$vr[,j], imaginary= z$vr[,j+1])
                vectors[,j+1] <- complex(real=z$vr[,j], imaginary=-z$vr[,j+1])
            }
        }
    }

    if( !only.values )
        return(list(values=values, vectors=vectors,alpha=alpha, beta=z$beta))
    else
        return(list(values=values, vectors=NULL, alpha=alpha, beta=z$beta))
}

# Generalized eigenvalue problem for symmetric pd matrix
# using Lapack dsygv
# A and B symmetric and B positive definite

geigen.dsygv <- function(A,B, only.values=FALSE) {

    # interface to Lapack dsygv
    # for generalized eigenvalue problem
    # symmetric (p.d.) matrices

    # if uplo=="U" upper triangular part of A, B contains upper triangular parts of matrix A, B
    # if uplo=="L" lower triangular part of A, B contains lower triangular parts of matrix A, B

    # ptype
    #   Specifies the problem type to be solved:
    #   = 1:  A*x = (lambda)*B*x
    #   = 2:  A*B*x = (lambda)*x
    #   = 3:  B*A*x = (lambda)*x

    ptype <- 1L
    # uplo <- match.arg(uplo)
    # uplo <- "L" #c("U","L")
    kuplo <- 2L

    dimA <- dim(A)
    n <- dimA[1]

    # jobev.char <- if(!only.values) "V" else "N"
    kjobev <- if(!only.values) 2L else 1L

    # calculate optimal workspace
    lwork <- -1L
    work  <- numeric(1)
    z <- .Fortran("xdsygv", as.integer(ptype), kjobev, kuplo, n, numeric(1), numeric(1),
                            numeric(1), work=work, lwork, info=integer(1L), PACKAGE="geigen")
    lwork <- as.integer(z$work[1])
    lwork <- max(lwork, as.integer(3*n-1))
    work  <- numeric(lwork)

    z <- .Fortran("xdsygv", as.integer(ptype), kjobev, kuplo, n, E=A, B,
                            w=numeric(n), work, lwork, info=integer(1L), PACKAGE="geigen")

    if( z$info != 0 ) .sygv_Lapackerror(z$info,n)

    if( !only.values )
        return(list(values=z$w, vectors=z$E, alpha=NULL, beta=NULL))
    else
        return(list(values=z$w, vectors=NULL, alpha=NULL, beta=NULL))
}

geigen.zggev <- function(A,B, only.values=FALSE) {

    # interface to Lapack zggev
    # for generalized eigenvalue problem
    # for symmetric matrices use zhegv

    # jobvl.char <- "N"
    # jobvr.char <- if(!only.values) "V" else "N"

    kjobvl <- 1L
    kjobvr <- if(!only.values) 2L else 1L

    dimA <- dim(A)
    n <- dimA[1]

    # calculate optimal workspace
    lwork <- -1L
    work  <- complex(1)
    rwork <- numeric(1)
    z <- .Fortran("xzggev", kjobvl, kjobvr, n, A, B, complex(1), complex(1),
                           complex(1), n, complex(1), n,
                           work=work, lwork, rwork, info=integer(1L), PACKAGE="geigen")
    lwork <- as.integer(Re(z$work[1]))
    lwork <- max(lwork, as.integer(2*n))
    work  <- complex(lwork)

    rwork <- numeric(8*n)
    tmp <- 0+0i

    if( !only.values )
        z <- .Fortran("xzggev", kjobvl, kjobvr, n, A, B, alpha=complex(n),
                               beta=complex(n), numeric(1),1L, vr=matrix(tmp,nrow=n,ncol=n), n,
                               work, lwork, rwork, info=integer(1L), PACKAGE="geigen")
    else
        z <- .Fortran("xzggev", kjobvl, kjobvr, n, A, B, alpha=complex(n),
                               beta=complex(n), numeric(1), 1L,  numeric(1), 1L,
                               work, lwork, rwork, info=integer(1L), PACKAGE="geigen")

    if( z$info != 0 ) .ggev_Lapackerror(z$info,n)

    # see comment in geigen.dggev
    values <- complexdiv(z$alpha,z$beta)

    if( !only.values )
        return(list(values=values, vectors=z$vr, alpha=z$alpha, beta=z$beta))
    else
        return(list(values=values, vectors=NULL, alpha=z$alpha, beta=z$beta))
}

# Generalized eigenvalue problem for hermitian matrix
# using Lapack zhegv
# A and B are Hermitian and B is positive definite

geigen.zhegv <- function(A,B, only.values=FALSE) {

    # A, B hermitian
    # B also positive definite
    # interface to Lapack zhegv
    # for generalized eigenvalue problem
    # symmetric (p.d.) matrices

    # if uplo=="U" upper triangular part of A, B contains upper triangular parts of matrix A, B
    # if uplo=="L" lower triangular part of A, B contains lower triangular parts of matrix A, B

    # ptype
    #   Specifies the problem type to be solved:
    #   = 1:  A*x = (lambda)*B*x
    #   = 2:  A*B*x = (lambda)*x
    #   = 3:  B*A*x = (lambda)*x

    ptype <- 1L
    # uplo <- match.arg(uplo)
    # uplo <- "L" #c("U","L")
    kuplo <- 2L

    dimA <- dim(A)
    n <- dimA[1]

    # jobev.char <- if(!only.values) "V" else "N"
    kjobev <- if(!only.values) 2L else 1L

    # calculate optimal workspace
    lwork <- -1L
    work  <- complex(1)
    z <- .Fortran("xzhegv", as.integer(ptype), kjobev, kuplo, n, complex(1), complex(1),
                           w=numeric(1), work=work, lwork,  numeric(1), info=integer(1L), PACKAGE="geigen")
    lwork <- as.integer(Re(z$work[1]))
    lwork <- max(lwork, as.integer(2*n-1))
    work  <- complex(lwork)

    rwork <- numeric(3*n-2)

    z <- .Fortran("xzhegv", as.integer(ptype), kjobev, kuplo, n, E=A, B,
                           w=numeric(n), work, lwork, rwork, info=integer(1L), PACKAGE="geigen")

    if( z$info != 0 ) .sygv_Lapackerror(z$info,n)

    if( !only.values )
        return(list(values=z$w, vectors=z$E, alpha=NULL, beta=NULL))
    else
        return(list(values=z$w, vectors=NULL, alpha=NULL, beta=NULL))
}
