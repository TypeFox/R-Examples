
#:::::::::::
#   dstup
#:::::::::::

subroutine  dstup (s, lds, nobs, nnull, qraux, jpvt, y, q, ldqr, ldqc, nq,_
                   info, work)

#  Purpose:  To perform QR decomposition of S=FR and to form F^{T}y, F^{T}QF's.

integer           lds, nobs, nnull, jpvt(*), ldqr, ldqc, nq, info
double precision  s(lds,*), y(*), qraux(*), q(ldqr,ldqc,*), work(*) 

#  On entry:
#      s          the S matrix spanning null space, of size (lds,nnull).
#      lds        leading dimension of s.
#      nobs       number of observations.
#      nnull      dimension of null space.
#      y          observations, of size (nobs).
#      q          the reproducing kernels, of size (ldqr,ldqc,nq).
#      ldqr       leading dimension for rows of q.
#      ldqc       leading dimension for columns of q.
#      nq         number of Q's.

#  On exit:
#      s,qraux,jpvt
#                 QR decomposition of S=FR.
#      y          F^{T} y.
#      q          F^{T}QF's.
#      info        0: normal termination.
#                 -1: dimension error.
#                 >0: rank(S)+1.

#  Work array:
#      work       of size at least (nobs).

#  Routines called directly:
#      Linpack -- dqrdc, dqrsl
#      Rkpack  -- dqrslm

#  Written:  Chong Gu, Statistics, Purdue, latest version 3/7/91.

double precision  dum
integer           j

info = 0

#   check dimension
if ( nobs < 1 | nobs > lds | nobs > ldqr | nobs > ldqc ) {
    info = -1
    return
}

#   QR decomposition of S=FR
    #   The indented line below is added on Mar 7, 1991, 
    #   with credit to Dick Franke
    for (j=1;j<=nnull;j=j+1)  jpvt(j) = 0
call  dqrdc (s, lds, nobs, nnull, qraux, jpvt, work, 1)

#   F^{T} y;  test rank of R
call  dqrsl (s, lds, nobs, nnull, qraux, y, dum, y, work, dum, dum, 01100,_
             info)
if ( info != 0 )  return

#   F^{T} Q_{*} F
for (j=1;j<=nq;j=j+1) {
    call  dqrslm (s, lds, nobs, nnull, qraux, q(1,1,j), ldqr, 0, info,_
                  work)
}

return
end

#...............................................................................

