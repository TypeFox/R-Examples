
# Provide as detailed as possible error messages for dgges/zgges/dggev/zggev/dsygv/zhegv

#    INFO is INTEGER
#      = 0:  successful exit
#      < 0:  if INFO = -i, the i-th argument had an illegal value.
#      = 1,...,N:
#            The QZ iteration failed.  (A,B) are not in Schur
#            form, but ALPHAR(j), ALPHAI(j), and BETA(j) should
#            be correct for j=INFO+1,...,N.
#      > N:  =N+1: other than QZ iteration failed in DHGEQZ.
#            =N+2: after reordering, roundoff changed values of
#                  some complex eigenvalues so that leading
#                  eigenvalues in the Generalized Schur form no
#                  longer satisfy SELCTG=.TRUE.  This could also
#                  be caused due to scaling.
#            =N+3: reordering failed in DTGSEN.

.gges_Lapackerror <- function(info, n) {

    if(info < 0 ) stop("(Internal) Bad call of Lapack routine")
    if(info > 0 && info <= n) {
        warning(paste0("QZ iteration failed but result should be correct for (alpha,beta) values[",info+1,":",n,"]"))
    } else if(info == n+1 ) {
        stop("Other than QZ iteration failed")
    } else if(info == n+2 ) {
        stop("Reordering inaccurate due to roundoff.")
    } else if(info == n+3 ) {
        stop("Complete failure of reordering")
    }
}

#    INFO is INTEGER
#      = 0:  successful exit
#      < 0:  if INFO = -i, the i-th argument had an illegal value.
#      = 1,...,N:
#            The QZ iteration failed.  No eigenvectors have been
#            calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)
#            should be correct for j=INFO+1,...,N.
#      > N:  =N+1: other than QZ iteration failed in DHGEQZ.
#            =N+2: error return from DTGEVC.

.ggev_Lapackerror <- function(info, n) {
    if(info < 0 ) stop("(Internal) Bad call of Lapack routine")
    if(info > 0 && info <= n) {
        warning(paste0("QZ iteration failed but result should be correct for eigenvalues[",info+1,":",n,"]"))
    } else if(info == n+1 ) {
        stop("Other than QZ iteration failed")
    } else if(info == n+2 ) {
        stop("Internal error.")
    }
}

#    INFO is INTEGER
#     = 0:  successful exit
#     < 0:  if INFO = -i, the i-th argument had an illegal value
#     > 0:  DPOTRF or DSYEV returned an error code:
#        <= N:  if INFO = i, DSYEV failed to converge;
#               i off-diagonal elements of an intermediate
#               tridiagonal form did not converge to zero;
#        > N:   if INFO = N + i, for 1 <= i <= N, then the leading
#               minor of order i of B is not positive definite.
#               The factorization of B could not be completed and
#               no eigenvalues or eigenvectors were computed.

.sygv_Lapackerror <- function(info, n) {
    if(info < 0 ) stop("(Internal) Bad call of Lapack routine")
    if(info > 0 && info <= n) {
        stop(paste("Convergence failure for", info," off-diagonal elements"))
    } else {
        stop(paste("Leading minor of order", info-n, "of B is not positive definite"))
    }
}

# .gges_Lapackerror(3,4)
# .ggev_Lapackerror(2,4)
# .sygv_Lapackerror(6,4)

.gsvd_Lapackerror <- function(info, n) {
    if(info < 0 ) stop("(Internal) Bad call of Lapack routine")
    if(info > 0 ) warning("Jacobi-type procedure failed to converge")
}
