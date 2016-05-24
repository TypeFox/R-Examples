
#:::::::::::
#   dcrdr
#:::::::::::

subroutine  dcrdr (s, lds, nobs, nnull, qraux, jpvt, q, ldq, nlaht,_
                   r, ldr, nr, cr, ldcr, dr, lddr, wk, info)

#  Purpose:  To compute auxiliary quantities cr and dr for posterior covariance

#  Usage:  Use s, qraux, jpvt, q, and nlaht returned by dsidr.

integer           lds, nobs, nnull, jpvt(*), ldq, ldr, nr, ldcr, lddr, info
double precision  s(lds,*), qraux(*), q(ldq,*), nlaht, r(ldr,*), cr(ldcr,*),_
                  dr(lddr,*), wk(2,*)

#  On entry:
#      s,qraux,jpvt
#                 QR-decomposition of  S = F R.
#      nobs       number of observations.
#      nnull      dimension of null space.
#      q          U^{T} F_{2}^{T} Q F_{2} U in BOTTOM-RIGHT corner's
#                     LOWER triangle and SUPER DIAGONAL;
#                 F_{2}^{T} Q F_{1} in BOTTOM-LEFT corner;
#      ldq        leading dimension of q.
#      nlaht      estimated log10(n*lambda).
#      r          R(t,s^{T}).
#      nr         length of s.

#  On exit:
#      cr         (M^{-1}-M^{-1}S(S^{T}M^{-1}S)^{-1}S^{T}M^{-1})R(t,s^{T})
#      dr         (S^{T}M^{-1}S)^{-1}S^{T}M^{-1}R(t,s^{T})
#      info        0: normal termination.
#                 >0: S is not of full rank: rank(S)+1 .
#                 -1: dimension error.
#                 -2: F_{2}^{T} Q F_{2} is not non-negative definite.
#      others     intact.

#  Work array:
#      wk         of size at least (2,nobs-nnull).

#  Routines called directly:
#      Blas    -- daxpy, dcopy, ddot
#      Linpack -- dpbfa, dpbsl, dqrsl, dtrsl
#      Other   -- dprmut, dset

#  Written:  Chong Gu, Statistics, Purdue, 2/27/96 at Ann Arbor.

double precision  dum, ddot
integer           i, j, n, n0

info = 0

#   check dimension
if ( nnull < 1 | nnull >= nobs | nobs > lds | nobs > ldq | ldr < nobs |_
     nr < 1 | ldcr < nobs | lddr < nnull ) {
    info = -1
    return
}

#   set working parameters
n0 = nnull
n = nobs - nnull

#   copy r to cr
for (j=1;j<=nr;j=j+1)  call  dcopy (nobs, r(1,j), 1, cr(1,j), 1)

#   compute  diag(I, U^{T}) F^{T} R(t,s^{T})
call  dcopy (n-2, q(n0+2,n0+1), ldq+1, wk, 1)
for (j=1;j<=nr;j=j+1) {
    call  dqrsl (s, lds, nobs, nnull, qraux, cr(1,j), dum, cr(1,j), dum,_
                 dum, dum, 01000, info)
    call  dqrsl (q(n0+2,n0+1), ldq, n-1, n-2, wk, cr(n0+2,j), dum,
                 cr(n0+2,j), dum, dum, dum, 01000, info)
}

#   compute  U ( T + n*lambdahat I )^{-1} diag(I, U^{T}) F^{T} R(t,s^{T})
call  dset (n, 10.d0 ** nlaht, wk(2,1), 2)
call  daxpy (n, 1.d0, q(n0+1,n0+1), ldq+1, wk(2,1), 2)
call  dcopy (n-1, q(n0+1,n0+2), ldq+1, wk(1,2), 2)
call  dpbfa (wk, 2, n, 1, info)
if ( info != 0 ) {
    info = -2
    return
}
for (j=1;j<=nr;j=j+1)  call  dpbsl (wk, 2, n, 1, cr(n0+1,j))
call  dcopy (n-2, q(n0+2,n0+1), ldq+1, wk, 1)
for (j=1;j<=nr;j=j+1)
    call  dqrsl (q(n0+2,n0+1), ldq, n-1, n-2, wk, cr(n0+2,j), cr(n0+2,j),_
                 dum, dum, dum, dum, 10000, info)

#   compute  dr
for (j=1;j<=nr;j=j+1) {
    for (i=1;i<=n0;i=i+1)
        dr(i,j) = cr(i,j) - ddot (n, cr(n0+1,j), 1, q(n0+1,i), 1)
    call  dtrsl (s, lds, n0, dr(1,j), 01, info)
    call  dprmut (dr(1,j), n0, jpvt, 1)
}

#   compute  cr
for (j=1;j<=nr;j=j+1) {
    call  dset (n0, 0.d0, cr(1,j), 1)
    call  dqrsl (s, lds, nobs, nnull, qraux, cr(1,j), cr(1,j),_
                 dum, dum, dum, dum, 10000, info)
}

return
end

#..............................................................................
