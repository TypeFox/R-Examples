
is.gsvd <- function(x) inherits(x,"xdgsvd") || inherits(x,"xzgsvd")

gsvd <- function(A,B) {

    if(!is.matrix(A)) stop("Argument A should be a matrix")
    if(!is.matrix(B)) stop("Argument B should be a matrix")

    dimA <- dim(A)
    dimB <- dim(B)

    if(dimA[1]==0) stop("Matrix A has zero rows/columns")
    if(dimB[1]==0) stop("Matrix B has zero rows/columns")
    if( dimA[2] != dimB[2] ) stop("Matrices A and B must have same number of columns")

    if(!all(is.finite(A))) stop("Matrix A may not contain infinite/NaN/NA")
    if(!all(is.finite(B))) stop("Matrix B may not contain infinite/NaN/NA")

    complex.A <- is.complex(A)
    complex.B <- is.complex(B)
    complex.AB <- FALSE
    if( is.complex(A) || is.complex(B) ) {
        # at least one matrix is complex ==> use complex routines
        if( !complex.A ) storage.mode(A) <- "complex"
        if( !complex.B ) storage.mode(B) <- "complex"
        complex.AB <- TRUE
    } else {
        # none of the matrices is complex ==> use double routines
        if( !is.double(A) ) storage.mode(A) <- "double"
        if( !is.double(B) ) storage.mode(B) <- "double"
    }

    kjobu <- 1L # Orthogonal matrix U is computed
    kjobv <- 1L # Orthogonal matrix V is computed
    kjobq <- 1L # Orthogonal matrix Q is computed

    m <- as.integer(dimA[1]) # number of rows A
    n <- as.integer(dimA[2]) # number of columns A and B
    p <- as.integer(dimB[1]) # number of rows B

    if( complex.AB ) {
        lwork <- -1L
        work  <- complex(1)
        z <- .Fortran("xzggsvd", kjobu, kjobv, kjobq,
                                 m, n, p,
                                 k=integer(1L), l=integer(1L),
                                 A=A, m, B=B, p,
                                 alpha=numeric(n), beta=numeric(n),
                                 U=matrix(complex(m*m),nrow=m), m, V=matrix(complex(p*p),nrow=p), p,
                                 Q=matrix(complex(n*n),nrow=n), n,
                                 work=work, lwork, rwork=numeric(2*n),
                                 iwork=integer(n),
                                 info=integer(1L), PACKAGE="geigen"
                      )

        lwork <- as.integer(z$work[1])
        work  <- complex(lwork)

        z <- .Fortran("xzggsvd", kjobu, kjobv, kjobq,
                                 m, n, p,
                                 k=integer(1L), l=integer(1L),
                                 A=A, m, B=B, p,
                                 alpha=numeric(n), beta=numeric(n),
                                 U=matrix(complex(m*m),nrow=m), m, V=matrix(complex(p*p),nrow=p), p,
                                 Q=matrix(complex(n*n),nrow=n), n,
                                 work=work, lwork, rwork=numeric(2*n),
                                 iwork=integer(n),
                                 info=integer(1L), PACKAGE="geigen"
                      )
        zclass <- "xzgsvd"
    } else {
        lwork <- -1L
        work  <- numeric(1)

        z <- .Fortran("xdggsvd", kjobu, kjobv, kjobq,
                                 m, n, p,
                                 k=integer(1L), l=integer(1L),
                                 A=A, m, B=B, p,
                                 alpha=numeric(n), beta=numeric(n),
                                 U=matrix(numeric(m*m),nrow=m), m, V=matrix(numeric(p*p),nrow=p), p,
                                 Q=matrix(numeric(n*n),nrow=n), n,
                                 work=work, lwork,
                                 iwork=integer(n),
                                 info=integer(1L), PACKAGE="geigen"
                      )

        lwork <- as.integer(z$work[1])
        work  <- numeric(lwork)

        z <- .Fortran("xdggsvd", kjobu, kjobv, kjobq,
                                 m, n, p,
                                 k=integer(1L), l=integer(1L),
                                 A=A, m, B=B, p,
                                 alpha=numeric(n), beta=numeric(n),
                                 U=matrix(numeric(m*m),nrow=m), m, V=matrix(numeric(p*p),nrow=p), p,
                                 Q=matrix(numeric(n*n),nrow=n), n,
                                 work=work, lwork,
                                 iwork=integer(n),
                                 info=integer(1L), PACKAGE="geigen"
                      )
        zclass <- "xdgsvd"
    }
    if(z$info!=0) .gsvd_Lapackerror(z$info)

    ret <- list(A=z$A, B=z$B, m=m, k=z$k, l=z$l, alpha=z$alpha, beta=z$beta, U=z$U, V=z$V, Q=z$Q)
    class(ret) <- zclass
    ret
}

# as it is coded in lapack's testing dgsvts.f
#makeR <- function(z) {
#    # generate R matrix for the Generalized SVD
#    # from dgsvts.f in Lapack folder testing/eig/dgstvs.f
#
#    k <- z$k
#    l <- z$l
#    m <- z$m
#    n <- ncol(z$A)
#    R <- matrix(0,nrow=k+l,ncol=k+l)
#    for( i in 1:min(k+l,m) ) {
#        for(j in i:(k+l) ) {
#            R[i,j] <- z$A[i,n-k-l+j]
#        }
#    }
#    if( m-k-l < 0 ) {
#        for(i in (m+1):(k+l) ) {
#            for( j in i:(k+l) ) {
#                R[i,j] <- z$B[i-k,n-k-l+j]
#            }
#        }
#    }
## the following two only if really required: not here
##    print(R)
##    if(n>k+l)  R <- cbind(matrix(0,nrow=(k+l),ncol=(n-k-l)),R)
#    R
#}

# See comments in Lapack source file dggsvd.f
# http://www.netlib.org/lapack/explore-html/dd/db4/dggsvd_8f.html
# http://www.netlib.org/lapack/lug/node36.html

# Extract the triangular matrix R from the gsvd object

gsvd.R <- function(z) {
    if( !is.gsvd(z) ) stop("gsvd.R only accepts an object created with function gsvd")

    m <- z$m
    k <- z$k
    l <- z$l
    n <- ncol(z$A)
    if( m - k - l < 0 ) {
       R <- (rbind(z$A,z$B[(m-k+1):l,])[1:(k+l),])
       if(n > k+l) R <- R[,-(1:(n-k-l))] # zero initial columns removed
    }
    else
        R <- z$A[1:(k+l), (n-k-l+1):n]
    R
}

# Extract the matrix [ 0 R ] from the gsvd object

gsvd.oR <- function(z) {
    if( !is.gsvd(z) ) stop("gsvd.oR only accepts an object created with function gsvd")

    k <- z$k
    l <- z$l
    n <- ncol(z$A)
    R <- gsvd.R(z)
    if( n > k+l ) oR <- cbind(matrix(0,nrow=k+l,ncol=n-k-l),R) else oR <- R
    oR
}

# Extract the matrix D1 from the gsvd object

gsvd.D1 <- function(z) {
    if( !is.gsvd(z) ) stop("gsvd.D1 only accepts an object created with function gsvd")

    m <- z$m
    k <- z$k
    l <- z$l

    D1 <- matrix(0,nrow=m, ncol=k+l)
    if(k > 0) D1[1:k,1:k] <- diag(k)
    if(m-k-l >= 0 ) {
        if(l > 0) D1[(k+1):(k+l),(k+1):(k+l)] <- diag(z$alpha[(k+1):(k+l)])
    } else {
        if(m > k) D1[(k+1):m,(k+1):m] <- diag(z$alpha[(k+1):m],nrow=m-k)
    }
    D1
}

# Extract the matrix D2 from the gsvd object

gsvd.D2 <- function(z) {
    if( !is.gsvd(z) ) stop("gsvd.D2 only accepts an object created with function gsvd")

    m <- z$m
    k <- z$k
    l <- z$l
    p <- nrow(z$B)
    D2 <- matrix(0,nrow=p, ncol=k+l)
    if(m-k-l >= 0 ) {
        if(l > 0) D2[1:l,(k+1):(k+l)] <- diag(z$beta[(k+1):(k+l)])
    } else { # m-k-l < 0 ===> k+l > m
        if(m > k) D2[1:(m-k), (k+1):m] <- diag(z$beta[(k+1):m],nrow=m-k)
        D2[(m-k+1):(l),(m+1):(k+l)] <- diag(k+l-m)
    }
    D2
}
