
#:::::::::::
#   dsytr
#:::::::::::

subroutine  dsytr (x, ldx, n, tol, info, work)

#  Acronym:  Double-precision SYmmetric matrix TRidiagonalization.

#  Purpose:  This routine performs the Householder tridiagonalization
#      algorithm on symmetric matrix `x', with truncation strategy as 
#      described in Gu, Bates, Chen, and Wahba (1988).
 
#  References:  1. Golud and Van Loan (1983) Matrix Computation. (pp.276-7)
#               2. Gu, Bates, Chen, and Wahba(1988), TR#823, Stat, UW-M.
#               3. Dongarra et al.(1979) LINPACK User's Guide. (Chap. 9)

#  Relation with LINPACK:  This routine computes the tridiagonalization
#      U^{T}XU=T, where X is symmetric, T is tridiagonal, and U is an 
#      orthogonal matrix as the product of Housholder matrices.  To compute 
#      U^{T}y or Uy for vector y, we can use routine `dqrsl' of LINPACK.  
#      The calling procedure is:
#
#        1. Create vector `qraux' by  
#             call  dcopy(n-2, x(2,1), ldx+1, qraux, 1)
#        2. Call `dqrsl' as
#             call  dqrsl (x(2,1), ldx, n-1, n-2, qraux, y(2), ... )
 
integer           ldx, n, info
double precision  x(ldx,*), tol, work(*)

#  On entry:
#      x          symmetric matrix, only LOWER triangle refered.
#      ldx        leading dimension of x.
#      n          size of matrix `x'.
#      tol        truncation tolarence; if zero, set to square machine
#                 precision.
 
#  On Exit:
#      x          diagonal:  diagonal elements of tridiag. transf.
#                 upper triangle:  off-diagonal of tridiag. transf.
#                 lower triangle:  overwritten by Householder factors.
#      info        0 :  normal termination.
#                 -1 :  dimension error.

#  Work array:
#      work       of size at least (n).

#  Routines called directly:
#      Fortran -- dfloat, dsqrt
#      Blas    -- daxpy, ddot, dscal
#      Blas2   -- dsymv, dsyr2

#  Written:  Chong Gu, Statistics, UW-Madison, latest version 8/29/88.

double precision  nrmtot, nrmxj, alph, toltot, tolcum, toluni, dn, ddot
integer           j

info = 0

#   check dimension
if ( ldx < n | n <= 2 ) {
    info = -1
    return
}

#   total Frobenius norm
nrmtot = ddot (n, x, ldx+1, x, ldx+1)
for ( j=1 ; j<n ; j=j+1 )  
    nrmtot = nrmtot + 2.d0 * ddot (n-j, x(j+1,j), 1, x(j+1,j), 1)

#   compute machine precision
toltot = 1.d0
while ( 1.d0 + toltot > 1.d0 )  toltot = toltot / 2.d0
toltot = 4.d0 * toltot ** 2

#   set truncation criterion
if ( toltot < tol )  toltot = tol
toltot = toltot * nrmtot
dn = dfloat (n)
toluni = toltot * 6.d0 / dn / ( dn - 1.d0 ) / ( 2.d0 * dn - 1.d0 )

#   initialization
tolcum = 0.d0

#   main process

for ( j=1 ; j<n-1 ; j=j+1 ) {
    #   deduct the F-norm of new diagonal element to update the remainder
    nrmtot = nrmtot - x(j,j) * x(j,j)

    #   compute norm of `b'
    nrmxj = ddot (n-j, x(j+1,j), 1, x(j+1,j), 1)

    #   cumulate the tolarence
    dn = dfloat (n-j)
    tolcum = tolcum + toluni * dn * dn

    #   set diagonal separation if truncation applicable
    if ( 2.d0 * nrmxj <= tolcum ) {     
        x(j,j+1) = 0.d0
        call  dscal (n-j, 0.d0, x(j+1,j), 1)
        #   deduct the norm truncated from the tolerance
        tolcum = tolcum - 2.d0 * nrmxj 
        toltot = toltot - 2.d0 * nrmxj
        next
    }

    #   Householder transform
    if ( x(j+1,j) < 0.d0 )  x(j,j+1) = dsqrt (nrmxj)
    else  x(j,j+1) = - dsqrt (nrmxj)
    nrmtot = nrmtot - 2.d0 * nrmxj

    #   b = sign(b_{1}) b / nrm(b) 
    call  dscal (n-j, -1.d0/x(j,j+1), x(j+1,j), 1)

    #   v = b + e_{1} 
    x(j+1,j) = 1.d0 + x(j+1,j)

    #   p = D v / v_{1} 
    alph = 1.d0 / x(j+1,j)
    call  dsymv ('l', n-j, alph, x(j+1,j+1), ldx, x(j+1,j), 1,_
                 0.d0, work(j+1), 1)

    #   w = p - (p^{T}v) v / (2 v_{1}) 
    alph = - ddot (n-j, work(j+1), 1, x(j+1,j), 1) / 2.d0 / x(j+1,j)
    call  daxpy (n-j, alph, x(j+1,j), 1, work(j+1), 1)

    #   D = D - v w^{T} - w v^{T} 
    call  dsyr2 ('l', n-j, -1.d0, x(j+1,j), 1, work(j+1), 1, x(j+1,j+1),_
                 ldx)

}

x(n-1,n) = x(n,n-1)

return
end

#...............................................................................

