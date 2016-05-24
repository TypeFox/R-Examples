
#:::::::::::
#   dgold
#:::::::::::

subroutine  dgold (vmu, q, ldq, n, z, low, upp, nlaht, score, varht,_
                   info, twk, work)

#  Purpose:  To evaluate GCV/GML function based on tridiagonal form and to
#      search minimum on an interval by golden section search.

character*1       vmu
integer           ldq, n, info
double precision  q(ldq,*), z(*), low, upp, nlaht, score, varht, twk(2,*),_
                  work(*)

#  On entry:
#      vmu        'v':  GCV criterion.
#                 'm':  GML criterion.
#                 'u':  unbiased risk estimate.
#      q          tidiagonal matrix in diagonal and super diagonal.
#      ldq        leading dimension of Q.
#      n          size of the matrix.
#      z          U^{T} F_{2}^{T} y.
#      low        lower limit of log10(n*lambda).
#      upp        upper limit of log10(n*lambda).
#      varht      known variance if vmu=='u'.

#  On exit:
#      nlaht      the estimated log(n*lambda).
#      score      the GCV/GML/URE score at the estimated lambda.
#      varht      the variance estimate at the estimated lambda.
#      info        0: normal termination.
#                 -1: dimension error.
#                 -2: tridiagonal form is not non-negative definite.
#                 -3: vmu is none of 'v', 'm', or 'u'.

#  Work arrays:
#      twk        of size at least (2,n).
#      work       of size at least (n).

#  Routines called directly:
#      Fortran -- dsqrt
#      Blas    -- daxpy, dcopy
#      Rkpack  -- dtrev
#      Other   -- dset

#  Written:  Chong Gu, Statistics, Purdue, latest version 12/29/91.

double precision  ratio, mlo, mup, tmpl, tmpu

ratio = ( dsqrt (5.d0) - 1.d0 ) / 2.d0

info = 0

#   interchange the boundaries if necessary
if ( upp < low ) {
    mlo = low
    low = upp
    upp = mlo
}

#   check vmu
if ( vmu != 'v' & vmu != 'm' & vmu != 'u' ) {
    info = -3
    return
}

#   check dimension
if ( n < 1 | n > ldq ) {
    info = -1
    return
}

#   initialize golden section search for scrht
mlo = upp - ratio * (upp - low)
call  dset (n, 10.d0 ** (mlo), twk(2,1), 2)
call  daxpy (n, 1.d0, q, ldq+1, twk(2,1), 2)
call  dcopy (n-1, q(1,2), ldq+1, twk(1,2), 2)
twk(1,1) = 10.d0**mlo
call  dtrev (vmu, twk, 2, n, z, tmpl, varht, info, work)
if ( info != 0 ) {
    info = -2
    return
}
mup = low + ratio * (upp - low)
call  dset (n, 10.d0 ** (mup), twk(2,1), 2)
call  daxpy (n, 1.d0, q, ldq+1, twk(2,1), 2)
call  dcopy (n-1, q(1,2), ldq+1, twk(1,2), 2)
twk(1,1) = 10.d0**mup
call  dtrev (vmu, twk, 2, n, z, tmpu, varht, info, work)
if ( info != 0 ) {
    info = -2
    return
}

#   golden section search for estimate of lambda
repeat {
    if ( mup - mlo < 1.d-7 )  break
    if ( tmpl < tmpu ) {
        upp = mup
        mup = mlo
        tmpu = tmpl
        mlo = upp - ratio * (upp - low)
        call  dset (n, 10.d0 ** (mlo), twk(2,1), 2)
        call  daxpy (n, 1.d0, q, ldq+1, twk(2,1), 2)
        call  dcopy (n-1, q(1,2), ldq+1, twk(1,2), 2)
        twk(1,1) = 10.d0**mlo
        call  dtrev (vmu, twk, 2, n, z, tmpl, varht, info, work)
        if ( info != 0 ) {
            info = -2
            return
        }
    }    
    else {
        low = mlo
        mlo = mup
        tmpl = tmpu
        mup = low + ratio * (upp - low)
        call  dset (n, 10.d0 ** (mup), twk(2,1), 2)
        call  daxpy (n, 1.d0, q, ldq+1, twk(2,1), 2)
        call  dcopy (n-1, q(1,2), ldq+1, twk(1,2), 2)
        twk(1,1) = 10.d0**mup
        call  dtrev (vmu, twk, 2, n, z, tmpu, varht, info, work)
        if ( info != 0 ) {
            info = -2
            return
        }
    }    
}

#   compute the return value
nlaht = ( mup + mlo ) / 2.d0
call  dset (n, 10.d0 ** (nlaht), twk(2,1), 2)
call  daxpy (n, 1.d0, q, ldq+1, twk(2,1), 2)
call  dcopy (n-1, q(1,2), ldq+1, twk(1,2), 2)
twk(1,1) = 10.d0**nlaht
call  dtrev (vmu, twk, 2, n, z, score, varht, info, work)
if ( info != 0 ) {
    info = -2
    return
}

return
end

#...............................................................................

