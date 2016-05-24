
#:::::::::::
#   dcore
#:::::::::::

subroutine  dcore (vmu, q, ldq, nobs, nnull, tol, z, job, limnla, nlaht,_
                   score, varht, info, twk, work)

#  Purpose:  To evaluate the GCV/GML score function at various trial values
#      of n*lambda using the tridiagonalization GCV/GML algorithm.  Perform
#      either golden section search or regular grid search for minimizing
#      n*lambda.

character*1       vmu
integer           ldq, nobs, nnull, job, info
double precision  q(ldq,*), tol, z(*), limnla(2), nlaht, score(*), varht,_
                  twk(2,*), work(*)

#  On entry:
#      vmu        'v':  GCV criterion.
#                 'm':  GML criterion.
#                 'u':  unbiased risk estimate.
#      q          F^{T} Q F, only refer the LOWER triangle of the BOTTOM-
#                 RIGHT corner, i.e., F_{2}^{T} Q F_{2}.
#      ldq        leading dimension of Q.
#      nobs       number of observations.
#      nnull      dimension of null space.
#      tol        tolerance of truncation.
#      z          F^{T} y.
#      job         0:  searching interval for nlaht chosen automatically.
#                 -1:  searching interval for nlaht provided by limnla.
#                 >0:  search regular grid points on [limnla(1),limnla(2)]:
#                        #(grids) = job + 1.
#      limnla     searching interval in log10 scale, see job.
#      varht      known variance if vmu=='u'.

#  On exit:
#      q          tridiagonal form in diagonal and superdiagonal of the
#                 corner, Householder factors in strict lower triangle of 
#                 the corner.
#      z          diag(I, U^{T}) F^{T} y.
#      limnla     see limnla of entry.
#      nlaht      the estimated log10(n*lambda).
#      score      job <= 0 :  the GCV/GML/URE score at nlaht.
#                 job  > 0 :  the GCV/GML/URE score at the regular grid points.
#      varht      variance estimate.
#      info        0 :  normal termination.
#                 -1 :  dimension error.
#                 -2 :  F_{2}^{T}QF_{2} is not non-negative definite.
#                 -3 :  vmu is none of 'v', 'm', or 'u'.

#  Work arrays:
#      twk        of size at least (2,nobs-nnull).
#      work       of size at least (nobs-nnull).

#  Routines called directly:
#      Fortran -- dfloat, dlog10
#      Blas    -- dasum, dcopy
#      Linpack -- dqrsl
#      Rkpack  -- deval, dgold, dsytr

#  Written:  Chong Gu, Statistics, Purdue, latest version 3/24/92.

double precision  dum, low, upp, dasum, mchpr
integer           n0, n, j

info = 0

#   check vmu
if ( vmu != 'v' & vmu != 'm' & vmu != 'u' ) {
    info = -3
    return
}

#   check dimension
if ( nnull < 1 | nobs <= nnull | nobs > ldq ) {
    info = -1
    return
}

#   set working parameters
n0 = nnull
n = nobs - nnull

#   tridiagonalization  U^{T} ( F_{2}^{T} Q F_{2} ) U = T
call  dsytr (q(n0+1,n0+1), ldq, n, tol, info, work)
             if ( info != 0 )  return

#   U^{T} z_{2}
call  dcopy (n-2, q(n0+2,n0+1), ldq+1, work, 1)
call  dqrsl (q(n0+2,n0+1), ldq, n-1, n-2, work, z(n0+2), dum, z(n0+2),_
             dum, dum, dum, 01000, info)

#   set searching range
if ( job == 0 ) {
    mchpr = 1.d0
    while ( 1.d0 + mchpr > 1.d0 )  mchpr = mchpr / 2.d0
    mchpr = mchpr * 2.d0
    limnla(2) = dmax1 (dasum (n, q(n0+1,n0+1), ldq+1) * 1.d2, mchpr)
    limnla(1) = limnla(2) * mchpr
    limnla(2) = dlog10 (limnla(2))
    limnla(1) = dlog10 (limnla(1))
}

low = limnla(1)
upp = limnla(2)
if ( job <= 0 ) {
    #   compute score and estimate nlaht thru golden-section search
    call dgold (vmu, q(n0+1,n0+1), ldq, n, z(n0+1), low, upp, nlaht,_
                score(1), varht, info, twk, work)
    if ( vmu == 'v' )  score(1) = score(1) * dfloat (nobs) / dfloat (n)
    if ( vmu == 'm' )  score(1) = score(1) * dfloat (n) / dfloat (nobs)
    if ( vmu == 'u' )  score(1) = score(1) * dfloat (n) / dfloat (nobs) + 2.d0 * varht
}
else {
    #   regular grid evaluation
    call  deval (vmu, q(n0+1,n0+1), ldq, n, z(n0+1), job, low, upp, nlaht,_
                 score, varht, info, twk, work)
    dum = dfloat (nobs) / dfloat (n)
    for (j=1;j<=job+1;j=j+1) {
        if ( vmu == 'v' )  score(j) = score(j) * dum
        if ( vmu == 'm' )  score(j) = score(j) / dum
        if ( vmu == 'u' )  score(j) = score(j) / dum + 2.d0 * varht
    }
}

return
end

#...............................................................................

