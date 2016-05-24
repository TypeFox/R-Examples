
#:::::::::::
#   dtrev
#:::::::::::

subroutine  dtrev (vmu, t, ldt, n, z, score, varht, info, work)

#  Acronym:  Double-precision TRidiagonal EValuation.

#  Purpose:  To compute the GCV/GML function and the related variance
#      estimate from the tridiagonal matrix `t' and data vector `z'.

#  References:  1. Gu, Bates, Chen, and Wahba(1988), TR#823, Stat, UW-M.
#               2. Dongarra et al. (1979) LINPACK User's Guide. (Chap. 4)

character*1       vmu
integer           n, info
double precision  t(ldt,*), z(*), score, varht, work(*)

#  On entry:
#      vmu        'v':  GCV.
#                 'm':  GML.
#                 'u':  unbiased risk estimate.
#      t          the positive definite tridiagonal matrix  T,
#                 stored in packed form:
#                     t(1,2:n):  off-diagonal
#                     t(2,1:n):  diagonal.
#      ldt        leading dimension of t.
#      n          the dimension of the matrix.
#      z          the appropriately transformed data vector.
#      varht      known variance if vmu=='u'.

#  On exit:
#      score      the GCV/GML/URE score.
#      varht      \hat\sigma^{2}.
#      info         -3:  vmu is none of 'v', 'm', or 'u'.
#                 > -3:  as from LINPACK's `dpbfa'.

#  Work array:
#      work       of size at least (n).

#  Routines called directly:
#      Fortran -- dexp, dfloat, dlog
#      Blas    -- dasum, dcopy, ddot, dscal
#      Linpack -- dpbfa, dpbsl

#  Written:  Chong Gu, Statistics, UW-Madison, latest version 12/29/91.

double precision  nume, deno, tmp, alph, la, dasum, ddot
integer           j

info = 0

#   check vmu
if ( vmu != 'v' & vmu != 'm' & vmu != 'u' ) {
    info = -3
    return
}

la = t(1,1)

#   standardize the matrix for numerical stability
alph = dfloat (n) / dasum (n, t(2,1), ldt)
call  dscal (n, alph, t(2,1), ldt)
call  dscal (n-1, alph, t(1,2), ldt)

#   decomposition
call  dpbfa (t, ldt, n, 1, info)
if ( info != 0 )  return

call  dcopy (n, z, 1, work, 1)
call  dpbsl (t, ldt, n, 1, work)
    
#   GCV computation
if ( vmu == 'v' ) {
    tmp = 1.d0 / t(2,n) / t(2,n)
    deno = tmp
    for (j=n-1;j>0;j=j-1) {
        tmp = ( 1.d0 + t(1,j+1) * t(1,j+1) * tmp ) / t(2,j) / t(2,j)
        deno = deno + tmp
    }
    nume = ddot (n, work, 1, work, 1) / dfloat (n)
    deno = deno / dfloat (n)
    varht = alph * la * nume / deno
    score = nume / deno / deno
}

#   GML computation
if ( vmu == 'm' ) {
    deno = dlog (t(2,n))
    for (j=n-1;j>0;j=j-1)  deno = deno + dlog (t(2,j))
    nume = ddot (n, z, 1, work, 1) / dfloat (n)
    varht = alph * la * nume
    score = nume * dexp (2.d0 * deno / dfloat (n))
}

#   unbiased risk computation
if ( vmu == 'u' ) {
    nume = ddot (n, work, 1, work, 1) / dfloat (n)
    tmp = 1.d0 / t(2,n) / t(2,n)
    deno = tmp
    for (j=n-1;j>0;j=j-1) {
        tmp = ( 1.d0 + t(1,j+1) * t(1,j+1) * tmp ) / t(2,j) / t(2,j)
        deno = deno + tmp
    }
    deno = deno / dfloat (n)
    score = alph * alph * la * la * nume - 2.d0 * varht * alph * la * deno
}

return
end

#...............................................................................
