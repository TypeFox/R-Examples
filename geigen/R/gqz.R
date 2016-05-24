

gqz <- function(A,B, sort=c("N","-","+","S","B","R")) {

    sort <- match.arg(sort)

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

    ksort <- match(sort, c("N","-","+","S","B","R"))
    # to be absolutely sure
    if(is.na(ksort)) stop("invalid sort argument")

    if(is.complex(A) || is.complex(B)) {
        if(!is.complex(A)) storage.mode(A) <- "complex"
        if(!is.complex(B)) storage.mode(B) <- "complex"
        geigen.xzgges(A,B, ksort)
    } else {
        # none of the matrices is complex ==> use double routines
        if( !is.double(A) ) storage.mode(A) <- "double"
        if( !is.double(B) ) storage.mode(B) <- "double"
        geigen.xdgges(A,B, ksort)
    }
}

geigen.xdgges <- function(A,B, ksort) {

    # interface to xdgges which calls Lapack dgges
    # for generalized eigenvalue problem
    # general real matrices

    # jobvsl.char <- "V"
    # jobvsr.char <- "V"

    kjobvsl <- 2L
    kjobvsr <- 2L

    dimA <- dim(A)
    n <- dimA[1]

    # calculate optimal workspace
    lwork <- -1L
    work  <- numeric(1)
    sdim  <- 0L
    bwork <- logical(n)

    z <- .Fortran("xdgges", kjobvsl, kjobvsr, ksort,
                            n, A, B, sdim, numeric(1), numeric(1),
                            numeric(1), numeric(1), numeric(1),
                            work=work, lwork, bwork, info=integer(1L), PACKAGE="geigen")
    lwork <- as.integer(z$work[1])
    lwork <- max(lwork, as.integer(8*n+16))
    work  <- numeric(lwork)

    z <- .Fortran("xdgges", kjobvsl, kjobvsr, ksort,
                            n, A=A, B=B, sdim=sdim, alphar=numeric(n), alphai=numeric(n),
                            beta=numeric(n), vsl=matrix(0,nrow=n,ncol=n), vsr=matrix(0,nrow=n,ncol=n),
                            work, lwork, bwork, info=integer(1L), PACKAGE="geigen")

    if( z$info != 0 ) .gges_Lapackerror(z$info,n)

    ret <- list(S=z$A, T=z$B, sdim=z$sdim, alphar=z$alphar, alphai=z$alphai, beta=z$beta, Q=z$vsl, Z=z$vsr)
    class(ret) <- "xdgges"
    ret
}

gevalues.xdgges <- function(x) {
    if(!inherits(x,"xdgges")) stop("gevalues only accepts an object constructed with gqz")
    
    if( all(x$alphai==0) ) {
        ret <- x$alphar/x$beta
    } else {
        ret <- complexdiv(complex(real=x$alphar,imaginary=x$alphai),x$beta)
    }
    ret
}

print.xdgges <- function(x, ...) {
    print(unclass(x), ...)
    invisible(x)
}

geigen.xzgges <- function(A,B, ksort) {

    # interface to xzgges which calls Lapack zgges
    # for generalized eigenvalue problem
    # general complex matrices  (A and B must be complex; tested above)

    # jobvsl.char <- "V"
    # jobvsr.char <- "V"

    kjobvsl <- 2L
    kjobvsr <- 2L

    dimA <- dim(A)
    n <- dimA[1]

    # calculate optimal workspace
    lwork <- -1L
    work  <- complex(1)
    sdim  <- 0L
    bwork <- logical(n)

    z <- .Fortran("xzgges", kjobvsl, kjobvsr, ksort,
                            n, A, B, sdim, complex(1), complex(1),
                            complex(1), complex(1),
                            work=work, lwork, numeric(1), bwork, info=integer(1L), PACKAGE="geigen")
    lwork <- as.integer(Re(z$work[1]))
    lwork <- max(lwork, as.integer(8*n+16))
    work  <- complex(lwork)

    rwork <- numeric(8*n)
    tmp <- 0+0i

    z <- .Fortran("xzgges", kjobvsl, kjobvsr, ksort,
                            n, A=A, B=B, sdim=sdim, alpha=complex(n),
                            beta=complex(n), vsl=matrix(tmp,nrow=n,ncol=n), vsr=matrix(tmp,nrow=n,ncol=n),
                            work, lwork, rwork, bwork, info=integer(1L), PACKAGE="geigen")

    if( z$info != 0 ) .gges_Lapackerror(z$info,n)

    ret <- list(S=z$A, T=z$B, sdim=z$sdim, alpha=z$alpha, beta=z$beta, Q=z$vsl, Z=z$vsr)
    class(ret) <- "xzgges"
    ret
}

print.xzgges <- function(x, ...) {
    print(unclass(x), ...) 
    invisible(x)
}

gevalues.xzgges <- function(x) {
    if(!inherits(x,"xzgges")) stop("gevalues only accepts an object constructed with gqz")
    
    ret <- complexdiv(x$alpha,x$beta)
}
