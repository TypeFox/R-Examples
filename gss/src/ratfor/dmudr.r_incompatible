#:::::::::::
#   dmudr
#:::::::::::

subroutine  dmudr (vmu,_
                   s, lds, nobs, nnull, q, ldqr, ldqc, nq, y,_     # inputs
                   tol, init, prec, maxite,_                       # tune para
                   theta, nlaht, score, varht, c, d,_              # outputs
                   wk, info)

#  Acronym:  Double precision MUltiple smoothing parameter DRiver.
 
#  Purpose:  This routine implements the iterative algorithm for minimizing
#      GCV/GML scores with multiple smoothing parameters described in 
#      Gu and Wahba (1991, SISSC).

#  WARNING:  Please be sure that you understand what this routine does before 
#      you call it.  Pilot runs with small problems are recommended.  This
#      routine performs VERY INTENSIVE numerical calculations for big nobs.

integer           lds, nobs, nnull, ldqr, ldqc, nq, init, maxite,_
                  info

double precision  s(lds,*), q(ldqr,ldqc,*), y(*), tol, prec,_
                  theta(*), nlaht, score, varht, c(*), d(*),_
                  wk(*)

character*1       vmu

#  On entry:
#      vmu        'v':  GCV criterion.
#                 'm':  GML criterion.
#                 'u':  unbiased risk estimate.
#      s          the matrix S, of size (lds,nnull).
#      nobs       the number of observations.
#      nnull      the dimension of the null space.
#      q          the matrices Q_{i}'s, of dimension (ldqr,ldqc,nq).
#      nq         the number of Q_{i}'s.
#      y          the response vector of size (nobs)
#      tol        the tolerance for truncation in the tridiagonalization; usually set to 0.d0.
#      init       0 :  no initial values provided for the theta.
#                 1 :  initial values provided for the theta.
#      theta      initial values of theta if init = 1.
#      prec       precision requested for the minimum score value, where precision
#                   is the weaker of the absolute and relative precisions.
#      maxite     maximum number of iterations allowed; usually 20 is enough.
#      varht      known variance if vmu=='u'.
 
#  On exit:
#      theta      the vector of parameter log10(theta) used in the final model,
#                 of dimension (nq).  -25 indicates effective minus infinity.
#      nlaht      the estimated  log10(n*lambda)|theta  in the final model.
#      score      the minimum GCV/GML/URE score found at (theta, nlaht).
#      varht      the variance estimate.
#      c,d        the coefficient estimates.
#      info        0 :  normal termination.
#                 -1 :  dimension error.
#                 -2 :  F_{2}^{T} Q_{*}^{theta} F_{2} !>= 0.
#                 -3 :  tuning parameters are out of scope.
#                 -4 :  fails to converge within maxite steps.
#                 -5 :  fails to find a reasonable descent direction.
#                 >0 :  the matrix S is rank deficient: rank(S)+1.
#      s,q,y      destroyed.
#      others     intact.

#  Work arrays:
#      wk         of size (nobs*nobs*(nq+2))

#  Routines called directly:
#      Rkpack  -- dmudr1

#  Routines called indirectly:
#      Blas    -- dasum, daxpy, dcopy, ddot, dnrm2, dscal, dswap, idamax
#      Blas2   -- dgemv, dsymv, dsyr2
#      Fortran -- dabs, dexp, dfloat, dlog, dlog10, dmax1, dsqrt
#      Linpack -- dpbfa, dpbsl, dpofa, dposl, dqrdc, dqrsl, dtrsl
#      Rkpack  -- dcoef, dcore, ddeev, deval, dgold, dmcdc, dqrslm,
#                 dstup, dsytr, dtrev
#      Other   -- dprmut, dset

#  Written:  Chong Gu, Statistics, Purdue, latest version 3/9/91.

integer  n, n0
integer  iqraux, itraux, itwk, iqwk, iywk, ithewk, ihes, igra, ihwk1, ihwk2,_
         igwk1, igwk2, ikwk, iwork1, iwork2, ijpvt, ipvtwk


n = nobs
n0 = nnull

iqraux = 1
itraux = iqraux + n0
itwk = itraux + (n-n0-2)
iqwk = itwk + 2 * (n-n0)
iywk = iqwk + n * n
ithewk = iywk + n
ihes = ithewk + nq
igra = ihes + nq * nq
ihwk1 = igra + nq
ihwk2 = ihwk1 + nq * nq
igwk1 = ihwk2 + nq * nq
igwk2 = igwk1 + nq
ikwk = igwk2 + nq
iwork1 = ikwk + (n-n0) * (n-n0) * nq
iwork2 = iwork1 + n
ijpvt = iwork2 + n
ipvtwk = ijpvt + n0

call  dmudr1 (vmu,_
             s, lds, nobs, nnull, q, ldqr, ldqc, nq, y,_     # inputs
             tol, init, prec, maxite,_                       # tune para
             theta, nlaht, score, varht, c, d,_              # outputs
             wk(iqraux), wk(ijpvt), wk(itwk), wk(itraux), wk(iqwk),_
             wk(iywk), wk(ithewk), wk(ihes), wk(igra), wk(ihwk1),_
             wk(ihwk2), wk(igwk1), wk(igwk2), wk(ipvtwk), wk(ikwk),_
             wk(iwork1), wk(iwork2),_
             info)

return
end
