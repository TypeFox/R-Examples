
#:::::::::::
#   deval
#:::::::::::

subroutine  deval (vmu, q, ldq, n, z, nint, low, upp, nlaht, score, varht,_
                   info, twk, work)

#  Purpose:  To evaluate GCV/GML function based on tridiagonal form and to
#      search minimum on an interval by equally spaced (in log10 scale) grid
#      search.

character*1       vmu
integer           ldq, n, nint, info
double precision  q(ldq,*), z(*), low, upp, nlaht, score(*), varht,_
                  twk(2,*), work(*)

#  On entry:
#      vmu        'v':  GCV criterion.
#                 'm':  GML criterion.
#                 'u':  unbiased risk estimate.
#      q          tidiagonal matrix in diagonal and super diagonal.
#      ldq        leading dimension of Q.
#      n          size of the matrix.
#      z          U^{T} F_{2}^{T} y.
#      nint       number of intervals (number of grids minus 1).
#      low        lower limit of log10(n*lambda).
#      upp        upper limit of log10(n*lambda).
#      varht      known variance if vmu=='u'.

#  On exit:
#      nlaht      the estimated log10(n*lambda).
#      score      the GCV/GML/URE score vector on grid points.
#      varht      the variance estimate at the estimated n*lambda.
#      info        0: normal termination.
#                 -1: dimension error.
#                 -2: tridiagonal form is not non-negative definite.
#                 -3: vmu or nint is out of scope.

#  Work arrays:
#      twk        array of length at least (2,n).
#      work       array of length at least (n).

#  Routines called directly:
#      Fortran -- dfloat
#      Blas    -- daxpy, dcopy
#      Rkpack  -- dtrev
#      Other   -- dset

#  Written:  Chong Gu, Statistics, Purdue, 12/29/91 latest version.

double precision  tmp, minscr, mlo, varhtwk
integer           j

info = 0

#   interchange boundaries if necessary
if ( upp < low ) {
    mlo = low
    low = upp
    upp = mlo
}

#   check job requests
if ( (vmu != 'v' & vmu != 'm' & vmu != 'u') | nint < 1 ) {
    info = -3
    return
}

#   check dimension
if ( 1 > n | n > ldq ) {
    info = -1
    return
}

#   evaluation
for (j=1;j<=nint+1;j=j+1) {
    tmp = low + dfloat (j-1) * ( upp - low ) / dfloat (nint)
    call  dset (n, 10.d0 ** (tmp), twk(2,1), 2)
    call  daxpy (n, 1.d0, q, ldq+1, twk(2,1), 2)
    call  dcopy (n-1, q(1,2), ldq+1, twk(1,2), 2)
    twk(1,1) = 10.d0**tmp
    call  dtrev (vmu, twk, 2, n, z, score(j), varht, info, work)
    if ( info != 0 ) {
        info = -2
        return
    }
    if ( score(j) <= minscr | j == 1 ) {
        minscr = score(j)
        nlaht = tmp
        varhtwk = varht
    }
}
varht = varhtwk

return
end

#...............................................................................

