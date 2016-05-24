
#::::::::::::
#   dsidr
#::::::::::::

subroutine  dsidr (vmu,_
                   s, lds, nobs, nnull, y, q, ldq,_       # data
                   tol, job, limnla,_                     # job requests
                   nlaht, score, varht, c, d,_            # output
                   qraux, jpvt, wk,_                      # work arrays
                   info)                                  # error message

#  Acronym:  Double precision SIngle smoothing parameter DRiver.
 
#  Purpose:  
#   
#      This routine is the double precision single smoothing parameter
#      driver of the RKPACK -- a minipackage for solving the equations
#      
#              ( n lambda I + Sigma ) c + S d  =  y
#                                        S' c  =  0  ,
#
#      where Sigma is n-by-n and S is n-by-M, and lambda is the so-called
#      smoothing parameter chosen to minimize the GCV criterion
#
#                            (1/n) || ( I - A(lambda) ) y || ** 2
#            V(lambda)  =   --------------------------------------   ,
#                            [ (1/n) tr ( I - A(lambda) ) ] ** 2
#
#      where A(lambda), satisfying
#
#                      A(lambda) y  =  Sigma c + S d  ,
#
#      is the so-called influence matrix, OR to minimize the GML criterion
#            
#                               (1/n) y' ( I - A(lambda) ) y
#            M(lambda)  =   ------------------------------------   ,
#                            det [ (I - A(lambda))+ ]^{1/(n-M)}
#
#      where det[(...)+] is the product of nonzero eigenvalues of (...).
#
#      The general theory behind this is described in Kimeldorf and Wahba
#      (1971), which seeks the minimizer of certain variational problem in 
#      reproducing kernel hilbert space.  The generalized cross validation
#      (GCV) method for choosing the smoothing parameter lambda is propos-
#      ed by Craven and Wahba (1979).  The GML criterion is described and
#      compared with the GCV by Wahba (1985).  An example of this general
#      scheme is the thin plate smoothing spline model, as described by
#      Wahba and Wendelberger (1980), and Bates et al. (1987).
#
#      RKPACK is the implementation of the GCV/GML algorithm based on the
#      Householder tridiagonalization, as proposed by Gu et al. (1988). 
#      It does not assume any structure of Sigma and S, except that S is
#      of full rank, Sigma is symmetric, and
#
#                  S' c  =  0   ===>   c' Sigma c  >=  0            (*)
#
#      The Sigma matrix is the reproducing kernel (or semi-kernel) evalu-
#      ated at the data points, and the matrix S is a set of null space 
#      basis evaluated at the data points.
#
#      Dsidr will do either golden-section search or regular grid search
#      for the minimizing lambda of V/M(lambda).  In the goden-section 
#      search case, it does assume bowl-shaped V/M(lambda) curve.  If this 
#      is not true, the user may specify shorter searching intervals on
#      which the curve may be bowl-shaped.  The precision of n*lambda is
#      1.d-7 in the log10 scale.  In the regular grid search case, it 
#      provides a "GCV/GML curve" on the searching interval.  (For the
#      later case user should provide `score' as a vector, though in the
#      golden section search case only minimum GCV/GML value is recorded.)
#
#      RKPACK is a cubic order package.  In fitting univariate smoothing
#      spline models, a linear order algorithm developed independently
#      by Hutchinson and deHoog (1985) and by O'Sullivan (1985) is recommended.
#      Code by Woltring (1986) and O'Sullivan is available from NETLIB.


character*1       vmu
integer           lds, nobs, nnull, ldq, job, jpvt(*), info
double precision  s(lds,*), y(*), q(ldq,*), tol, limnla(2), nlaht, score(*),_
                  varht, c(*), d(*), qraux(*), wk(*)


#  On entry:
#      vmu        'v':  GCV criterion.
#                 'm':  GML criterion.
#                 'u':  unbiased risk estimate.
#      s          the matrix S of size (nobs,nnull).
#      lds        the leading dimension of s.
#      nobs       the number of observations.
#      nnull      the dimension of the null space.
#      y          the observations.
#      q          the matrix Q, only the lower triangle referred.
#      tol        tolerance for truncation in `dsytr'.  If 0.d0, set to 
#                 square of machine precision.
#      job        <=0 : golden-section search
#                     0 --  searching interval specified automatically.
#                    -1 --  search on (limnla(1), limnla(2)).
#                  >0 : regular grid search on [limnla(1), limnla(2)]
#                     #(grids) = job + 1.
#      limnla     the searching interval (in log10 scale), see job.
#      varht      known variance if vmu=='u'.

#  On exit:
#      nlaht      the GCV/GML/URE estimate of log10(nobs*lambda).
#      limnla     searching range for nlaht.
#      score      job <= 0 :  GCV/GML/URE value at nlaht.
#                 job >  0 :  GCV/GML/URE vector on the regular grid points.
#      varht      the variance estimate.
#      c          the parameters c.
#      d          the parameters d.
#      s,qraux,jpvt
#                 QR decomposition of S=FR, as from Linpack `dqrdc'.
#      q          first nnull columns: F^{T} Q F_{1}.
#                 BOTTOM-RIGHT corner: tridiagonalization of 
#                                      F_{2}^{T} Q F_{2}.
#      info        0: normal termination.
#                 -1: dimension error.
#                 -2: F_{2}^{T} Q F_{2} !>= 0.
#                 -3: vmu is out of scope.
#                 >0: the matrix S is rank deficient: rank(S)+1.
#      others     intact.

#  Work arrays:
#      wk         of size at least (3*nobs).

#  Routines called directly:
#      Rkpack  -- dcoef, dcore, dstup

#  Routines called indirectly:
#      Fortran -- dexp, dfloat, dlog, dlog10, dsqrt
#      Blas    -- dasum, daxpy, dcopy, ddot, dnrm2, dscal
#      Blas2   -- dsymv, dsyr2
#      Linpack -- dpbfa, dpbsl, dqrdc, dqrsl, dtrsl
#      Rkpack  -- deval, dgold, dqrslm, dsytr, dtrev
#      Other   -- dprmut, dset

#  Written:  Chong Gu, Statistics, Purdue, latest version 12/29/91.


info = 0

#   check dimension
if ( nnull < 1 | nnull >= nobs | nobs > lds | nobs > ldq ) {
    info = -1
    return
}

#   check vmu
if ( vmu != 'v' & vmu != 'm' & vmu != 'u' ) {
    info = -3
    return
}

#   main process

call  dstup (s, lds, nobs, nnull, qraux, jpvt, y, q, ldq, nobs, 1, info,_
             wk)
if ( info != 0 )  return

call  dcore (vmu, q, ldq, nobs, nnull, tol, y, job, limnla, nlaht, score,_
             varht, info, wk, wk(2*nobs+1))
if ( info != 0 )  return

call  dcoef (s, lds, nobs, nnull, qraux, jpvt, y, q, ldq, nlaht, c, d,_
             info, wk)

return
end

#...............................................................................
