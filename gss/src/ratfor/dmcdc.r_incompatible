
#:::::::::::
#   dmcdc
#:::::::::::

subroutine  dmcdc (a, lda, p, e, jpvt, info)

#  Acronym:  Double precision Modified Cholesky DeComposition.

#  Purpose:  This routine implements the modified Cholesky decomposition as
#      described by Gill, Murray, and Wright (p.111, Practical Optimization,
#      Academic Press, 1981).  The parameter delta is set to the maximum of
#      1.d-7 * (average diag) and 1.d-10.  Pivoting is enforced.  The result
#      is compatible with the Linpack routine `dposl'.

integer           lda, p, jpvt(*), info
double precision  a(lda,*), e(*)

#  On entry:
#      a          a symmetric matrix in the UPPER triangle.
#      lda        the leading dimension of a.
#      p          the size of a.

#  On exit:
#      a          the Cholesky factor  R  of  P^{T}AP + E = R^{T} R, where P 
#                 is a permutation matrix.
#      e          the amount of diagonal modification, of size (p).
#      jpvt       the permutation P, jpvt(j) contains the index of diagonal 
#                 element moved to j-th position, of size (p).
#      info        0:  normal termination.
#                 -1:  dimension error.

#  Routines called directly:
#      Blas    -- dasum, ddot, dscal, dswap, idamax
#      Fortran -- dabs, dmax1, dfloat, dsqrt

#   Written:  Chong Gu, Statistics, UW-Madison, latest version 9/16/88.

double precision  beta, delta, theta, tmp, dasum, ddot
integer           i, j, jmax, jtmp, idamax

info = 0

#   check dimension
if ( lda < p | p < 1 ) {
    info = -1
    return
}

#   compute constants
tmp = 1.d0
while ( 1.d0 + tmp > 1.d0 )  tmp = tmp / 2.d0
jmax = idamax (p, a, lda+1)
beta = dmax1 (2.d0 * tmp, dabs (a(jmax,jmax)))
tmp = dsqrt (dfloat (p*p-1))
if ( tmp < 1.d0 )  tmp = 1.d0
for (j=2;j<=p;j=j+1) {
    jmax = idamax (j-1, a(1,j), 1)
    beta = dmax1 (beta, dabs (a(jmax,j)) / tmp)
}
delta = dasum (p, a, lda+1) / dfloat (p) * 1.d-7
delta = dmax1 (delta, 1.d-10)
for (j=1;j<=p;j=j+1)  jpvt(j) = j

#   compute  P^{T}AP + E = LDL^{T}
for (j=1;j<=p;j=j+1) {
    #   pivoting
    jmax = idamax (p-j+1, a(j,j), lda+1) + j - 1
    if ( jmax != j ) {
        call  dswap (j-1, a(1,j), 1, a(1,jmax), 1)
        call  dswap (jmax-j-1, a(j,j+1), lda, a(j+1,jmax), 1)
        call  dswap (p-jmax, a(j,jmax+1), lda, a(jmax,jmax+1), lda)
        tmp = a(j,j)
        a(j,j) = a(jmax,jmax)
        a(jmax,jmax) = tmp
        jtmp = jpvt(j)
        jpvt(j) = jpvt(jmax)
        jpvt(jmax) = jtmp
    }
    #   compute j-th column of L^{T}
    for (i=1;i<j;i=i+1)  a(i,j) = a(i,j) / a(i,i)
    #   update j-th row and determine the parameter theta
    for (i=j+1;i<=p;i=i+1)
        a(j,i) = a(j,i) - ddot (j-1, a(1,j), 1, a(1,i), 1)
    #   specify theta
    if ( j == p )  theta = 0.d0
    else {
        jmax = idamax (p-j, a(j,j+1), lda) + j
        theta = dabs (a(j,jmax))
    }
    #   compute diagonal d(j,j)
    tmp = dmax1 (delta, dabs (a(j,j)), theta ** 2 / beta)
    e(j) = tmp - a(j,j)
    a(j,j) = tmp
    #   update remaining diagonals
    for (i=j+1;i<=p;i=i+1)  a(i,i) = a(i,i) - a(j,i) ** 2 / a(j,j)
}

#   scale
for (j=1;j<=p;j=j+1) {
    a(j,j) = dsqrt (a(j,j))
    call  dscal (p-j, a(j,j), a(j,j+1), lda)
}

return
end

#...............................................................................
