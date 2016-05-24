
#:::::::::::
#   dcoef
#:::::::::::

subroutine  dcoef (s, lds, nobs, nnull, qraux, jpvt, z, q, ldq, nlaht,_
                   c, d, info, twk)

#  Purpose:  To compute the estimated coefficients of the model.

integer           lds, nobs, nnull, jpvt(*), ldq, info
double precision  s(lds,*), qraux(*), z(*), q(ldq,*), nlaht, c(*), d(*),_
                  twk(2,*)

#  On entry:
#      s,qraux,jpvt
#                 QR-decomposition of  S = F R.
#      lds        leading dimension of s.
#      nobs       number of observations.
#      nnull      dimension of null space.
#      z          diag(I, U^{T}) F^{T} y.
#      q          U^{T} F_{2}^{T} Q F_{2} U in BOTTOM-RIGHT corner's
#                     LOWER triangle and SUPER DIAGONAL;
#                 F_{2}^{T} Q F_{1} in BOTTOM-LEFT corner.
#      ldq        leading dimension of q.
#      nlaht      estimated log10(n*lambda).

#  On exit:
#      c          parameters c.
#      d          parameters d.
#      info        0: normal termination.
#                 >0: S is not of full rank: rank(S)+1 .
#                 -1: dimension error.
#                 -2: F_{2}^{T} Q F_{2} is not non-negative definite.

#  Work array:
#      twk        of size at least (2,nobs-nnull).

#  Routines called directly:
#      Blas    -- daxpy, dcopy, ddot
#      Linpack -- dpbfa, dpbsl, dqrsl, dtrsl
#      Other   -- dprmut, dset

#  Written:  Chong Gu, Statistics, UW-Madison, 5/4/88 at Yale.

double precision  dum, ddot
integer           n, n0

info = 0

#   check dimension
if ( nnull < 1 | nnull >= nobs | nobs > lds | nobs > ldq ) {
    info = -1
    return
}

#   set working parameters
n0 = nnull
n = nobs - nnull

#   compute  U ( T + n*lambdahat I )^{-1} z
call  dset (n, 10.d0 ** nlaht, twk(2,1), 2)
call  daxpy (n, 1.d0, q(n0+1,n0+1), ldq+1, twk(2,1), 2)
call  dcopy (n-1, q(n0+1,n0+2), ldq+1, twk(1,2), 2)
call  dpbfa (twk, 2, n, 1, info)
if ( info != 0 ) {
    info = -2
    return
}
call  dpbsl (twk, 2, n, 1, z(n0+1))
call  dcopy (n-2, q(n0+2,n0+1), ldq+1, twk, 1)
call  dqrsl (q(n0+2,n0+1), ldq, n-1, n-2, twk, z(n0+2), z(n0+2), dum,_
             dum, dum, dum, 10000, info)

#   compute  c
call  dset (n0, 0.d0, c, 1)
call  dcopy (n, z(n0+1), 1, c(n0+1), 1)
call  dqrsl (s, lds, nobs, nnull, qraux, c, c, dum, dum, dum, dum, 10000,_
             info)

#   compute  d
for (j=1;j<=n0;j=j+1) {
    d(j) = z(j) - ddot (n, z(n0+1), 1, q(n0+1,j), 1)
}
call  dtrsl (s, lds, n0, d, 01, info)
call  dprmut (d, n0, jpvt, 1)

return
end

#...............................................................................

