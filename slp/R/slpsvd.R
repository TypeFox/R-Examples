.slpsvd <- function(A, N, K) {

    ## If N is passed in as floating point, the cast to 
    ## as.integer() in the Fortran call does not quite work properly, so
    ## force it to integer now, then check it against the square matrix A
    if(!is.integer(N)) {
      N<-as.integer(floor(N));
    } 
    stopifnot(dim(A)[1] == dim(A)[2], dim(A)[1] == N)

    out <- .Fortran("SLPSVD", M = as.integer(N), N = as.integer(N), 
                    A = as.double(A), LDA = as.integer(N), S = double(N), 
                    VT = double(N * N), LDVT = as.integer(N), 
                    WORK = double(3 * N * N + 7 * N), LWORK = as.integer(3 * N * N + 7 * N), 
                    IWORK = integer(8 * N), NLSV = as.integer(K), 
                    PACKAGE='slp')

    out <- list(d = out$S, u = matrix(out$A, nrow = N, ncol = N)[, 1:K])
    return(out)
}

