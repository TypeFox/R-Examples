 
#:::::::::::
#   dsms
#:::::::::::

subroutine  dsms (s, lds, nobs, nnull, jpvt, q, ldq, nlaht,_
                  sms, ldsms, wk, info)

#  Purpose:  To compute the auxiliary quantity sms for posterior covariance

#  Usage:  Use s, qraux, jpvt, q, and nlaht returned by dsidr.

integer           lds, nobs, nnull, jpvt(*), ldq, ldsms, info
double precision  s(lds,*), q(ldq,*), nlaht, sms(ldsms,*), wk(2,*)

#  On entry:
#      s,jpvt     QR-decomposition of  S = F R.
#      nobs       number of observations.
#      nnull      dimension of null space.
#      q          U^{T} F_{2}^{T} Q F_{2} U in BOTTOM-RIGHT corner's
#                     LOWER triangle and SUPER DIAGONAL;
#                 F_{2}^{T} Q F_{1} in BOTTOM-LEFT corner;
#                 F_{1}^{T} Q F_{1} in UPPER-LEFT corner's LOWER triangle.
#      ldq        leading dimension of q.
#      nlaht      estimated log10(n*lambda).

#  On exit:
#      sms        (S^{T}M^{-1}S)^{-1}.
#      info        0: normal termination.
#                 >0: S is not of full rank: rank(S)+1 .
#                 -1: dimension error.
#                 -2: F_{2}^{T} Q F_{2} is not non-negative definite.
#      inputs     intact but UPPER-RIGHT corner of q was used as work array.

#  Work array:
#      wk         of size at least (2,nobs-nnull).

#  Routines called directly:
#      Blas    -- daxpy, dcopy, ddot
#      Linpack -- dpbfa, dpbsl, dqrsl, dtrsl
#      Other   -- dprmut, dset

#  Written:  Chong Gu, Statistics, Purdue, latest version 4/17/92.

double precision  dum, ddot
integer           i, j, n, n0

info = 0

#   check dimension
if ( nnull < 1 | nnull >= nobs | nobs > lds | nobs > ldq | ldsms < nnull ) {
    info = -1
    return
}

#   set working parameters
n0 = nnull
n = nobs - nnull

#   compute  sms

#   U^{T} (F_{2}^{T} Q F_{1})
call  dcopy (n-2, q(n0+2,n0+1), ldq+1, wk, 1)
for (j=1;j<=n0;j=j+1) {
    call  dcopy (n, q(n0+1,j), 1, q(j,n0+1), ldq)
    call  dqrsl (q(n0+2,n0+1), ldq, n-1, n-2, wk, q(n0+2,j), dum,
                 q(n0+2,j), dum, dum, dum, 01000, info)
}
#   U^{T} (F_{2}^{T}QF_{2} + n lambda I)^{-1} (F_{2}^{T}QF_{1})
call  dset (n, 10.d0 ** nlaht, wk(2,1), 2)
call  daxpy (n, 1.d0, q(n0+1,n0+1), ldq+1, wk(2,1), 2)
call  dcopy (n-1, q(n0+1,n0+2), ldq+1, wk(1,2), 2)
call  dpbfa (wk, 2, n, 1, info)
if ( info != 0 ) {
    info = -2
    return
}
for (j=1;j<=n0;j=j+1)  call  dpbsl (wk, 2, n, 1, q(n0+1,j))
#   (F_{2}^{T}QF_{2} + n lambda I)^{-1} (F_{2}^{T}QF_{1})
call  dcopy (n-2, q(n0+2,n0+1), ldq+1, wk, 1)
for (j=1;j<=n0;j=j+1) {
    call  dqrsl (q(n0+2,n0+1), ldq, n-1, n-2, wk, q(n0+2,j), q(n0+2,j),_
                 dum, dum, dum, dum, 10000, info)
}
#   (F_{1}^{T}QF_{1} + n lambda I) -
#   (F_{1}^{T}QF_{2}^{T}) (F_{2}^{T}QF_{2} + n lambda I)^{-1} (F_{2}^{T}QF_{1})
for (i=1;i<=n0;i=i+1) {
    for (j=1;j<i;j=j+1)  sms(i,j) = sms(j,i)
    for (j=i;j<=n0;j=j+1)
        sms(i,j) = q(j,i) - ddot (n, q(n0+1,j), 1, q(i,n0+1), ldq)
    sms(i,i) = sms(i,i) + 10.d0**nlaht
}
#   R^{-1} ... R^{-T} and permutation
for (j=1;j<=n0;j=j+1)  call  dtrsl (s, lds, n0, sms(1,j), 01, info)
for (i=1;i<=n0;i=i+1) {
    call  dcopy (n0, sms(i,1), ldsms, wk, 1)
    call  dtrsl (s, lds, n0, wk, 01, info)
    call  dprmut (wk, n0, jpvt, 1)
    call  dcopy (n0, wk, 1, sms(i,1), ldsms)
}
for (j=1;j<=n0;j=j+1)  call  dprmut (sms(1,j), n0, jpvt, 1)

#   restore  F_{2}^{T} Q F_{1} to the BOTTOM-LEFT corner of q
for (j=1;j<=n0;j=j+1)  call  dcopy (n, q(j,n0+1), ldq, q(n0+1,j), 1)

return
end

#..............................................................................
