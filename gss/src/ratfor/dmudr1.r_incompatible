#:::::::::::
#   dmudr1
#:::::::::::

subroutine  dmudr1 (vmu,_
                   s, lds, nobs, nnull, q, ldqr, ldqc, nq, y,_     # inputs
                   tol, init, prec, maxite,_                       # tune para
                   theta, nlaht, score, varht, c, d,_              # outputs
                   qraux, jpvt, twk, traux, qwk, ywk, thewk,_      # work arrays
                   hes, gra, hwk1, hwk2, gwk1, gwk2, pvtwk,_
                   kwk, work1, work2,_
                   info)

#  Acronym:  Double precision MUltiple smoothing parameter DRiver.
 
#  Purpose:  This routine implements the iterative algorithm for minimizing
#      GCV/GML scores with multiple smoothing parameters described in 
#      Gu and Wahba(1988, Minimizing GCV/GML scores with multiple smoothing
#      parameters via the Newton method).

#  WARNING:  Please be sure that you understand what this routine does before 
#      you call it.  Pilot runs with small problems are recommended.  This
#      routine performs VERY INTENSIVE numerical calculations for big nobs.

integer           lds, nobs, nnull, ldqr, ldqc, nq, init, maxite,_
                  jpvt(*), pvtwk(*), info

double precision  s(lds,*), q(ldqr,ldqc,*), y(*), tol, prec,_
                  theta(*), nlaht, score, varht, c(*), d(*),_
                  qraux(*), traux(*), twk(2,*), qwk(ldqr,*), ywk(*),_
                  thewk(*), hes(nq,*), gra(*), hwk1(nq,*), hwk2(nq,*),_
                  gwk1(*), gwk2(*), kwk(nobs-nnull,nobs-nnull,*),_
                  work1(*), work2(*)

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
#      tol        the tolerance for truncation in the tridiagonalization.
#      init       0 :  No initial values provided for the theta.
#                 1 :  Initial values provided for the theta.
#      theta      initial values of theta if init = 1.
#      prec       precision requested for the minimum score value.
#      maxite     maximum number of iterations allowed.
#      varht      known variance if vmu=='u'.
 
#  On exit:
#      theta      the vector of parameter log10(theta) used in the final model,
#                 of dimension (nq).  -25 indicates effective minus infinity.
#      nlaht      the estimated  log10(n*lambda)|theta  in the final model.
#      score      the minimum GCV/GML/URE score found at (theta, nlaht).
#      varht      the variance estimate.
#      c,d        the coefficients estimates.
#      info        0 :  normal termination.
#                 -1 :  dimension error.
#                 -2 :  F_{2}^{T} Q_{*}^{theta} F_{2} !>= 0.
#                 -3 :  tuning parameters are out of scope.
#                 -4 :  fails to converge within maxite steps.
#                 -5 :  fails to find a reasonable descent direction.
#                 >0 :  the matrix S is rank deficient: rank(S)+1.

#  Work arrays:
#      qraux      of size at least (nnull).
#      jpvt       of size at least (nnull).
#      twk        of size at least (2,nobs-nnull).
#      traux      of size at least (nobs-nnull-2).
#      qwk        of size at least (nobs,nobs).
#      ywk        of size at least (nobs).
#      thewk      of size at least (nq).
#      hes        of size at least (nq,nq).
#      gra        of size at least (nq).
#      hwk1-2     of sizes at least (nq,nq).
#      gwk1-2     of sizes at least (nq).
#      pvtwk      of size at least (nq).
#      kwk        of size at least (nobs-nnull,nobs-nnull,nq).
#      work1-2    of sizes at least (nobs).

#  Routines called directly:
#      Blas    -- dasum, daxpy, dcopy, ddot, dscal, idamax
#      Blas2   -- dsymv
#      Fortran -- dfloat, dlog, dlog10
#      Linpack -- dpofa, dposl, sqrsl
#      Rkpack  -- dcoef, dcore, ddeev, dmcdc, dstup
#      Other   -- dprmut, dset

#  Routines called indirectly:
#      Blas    -- dasum, daxpy, dcopy, ddot, dnrm2, dscal, dswap, idamax
#      Blas2   -- dgemv, dsymv, dsyr2
#      Fortran -- dabs, dexp, dfloat, dlog, dlog10, dsqrt
#      Linpack -- dpbfa, dpbsl, dqrdc, dqrsl, dtrsl
#      Rkpack  -- deval, dgold, dqrslm, dsytr, dtrev
#      Other   -- dprmut, dset

#  Written:  Chong Gu, Statistics, Purdue, latest version 1/6/92.

double precision  alph, scrold, scrwk, nlawk, limnla(2),_
                  tmp, dasum, ddot
integer           n, n0, i, j, iwk, maxitwk, idamax, job

info = 0

#   set working parameters
n0 = nnull
n = nobs - nnull
maxitwk = maxite

#   check tuning parameters
if ( (vmu != 'v' & vmu != 'm' & vmu != 'u') | (init != 0 & init != 1) |_
     (maxitwk <=0) | (prec <= 0.d0) ) {
    info = -3
    return
}

#   check dimension
if ( lds < nobs | nobs <= n0 | n0 < 1 | ldqr < nobs | ldqc < nobs |_
     nq <= 0 ) {
    info = -1
    return
}

#   initialize
call  dstup (s, lds, nobs, n0, qraux, jpvt, y, q, ldqr, ldqc, nq, info,_
             work1)
if ( info != 0 )  return
if ( init == 1 )  call  dcopy (nq, theta, 1, thewk, 1)
else {
#   use the "plug-in" weights as the starting theta
    for (i=1;i<=nq;i=i+1) {
        thewk(i) = dasum (n, q(n0+1,n0+1,i), ldqr+1)
        if ( thewk(i) > 0.d0 )  thewk(i) = 1.d0 / thewk(i)
    }
    #   fit an initial model
    for (j=1;j<=nobs;j=j+1)  call  dset (nobs-j+1, 0.d0, qwk(j,j), 1)
    for (i=1;i<=nq;i=i+1) {
        for (j=1;j<=nobs;j=j+1)
            call  daxpy (nobs-j+1, thewk(i), q(j,j,i), 1, qwk(j,j), 1)
    }
    call  dcopy (nobs, y, 1, ywk, 1)
    call  dcore (vmu, qwk, ldqr, nobs, n0, tol, ywk, 0, limnla, nlawk,_
                 scrwk, varht, info, twk, work1)
    if (info != 0 )  return
    call  dcoef (s, lds, nobs, n0, qraux, jpvt, ywk, qwk, ldqr, nlawk,_
                 c, d, info, twk)
    #   assign weights due to norm  \theta^{2}c^{T}(Q_{i})c
    call  dqrsl (s, lds, nobs, n0, qraux, c, tmp, c, tmp, tmp, tmp,_
                 01000, info)
    for (i=1;i<=nq;i=i+1) {
        call  dsymv('l', n, thewk(i), q(n0+1,n0+1,i), ldqr, c(n0+1), 1,_
                    0.d0, work1, 1)
        thewk(i) = ddot (n, c(n0+1), 1, work1, 1) * thewk(i)
        if ( thewk(i) > 0.d0 )  thewk(i) = dlog10 (thewk(i))
        else  thewk(i) = -25.d0
    }
}
scrold = 1.d10

#   main process

job = 0
repeat {
    #   nq == 1
    if ( nq == 1 ) {
        theta(1) = 0.d0
        break
    }
    #   form  Qwk = \sum_{i=1}^{nq} \thewk_{i} Q_{i}
    for (j=1;j<=nobs;j=j+1)  call  dset (nobs-j+1, 0.d0, qwk(j,j), 1)
    for (i=1;i<=nq;i=i+1) {
        if ( thewk(i) <= -25.d0 )  next
        for (j=1;j<=nobs;j=j+1)
            call  daxpy (nobs-j+1, 10.d0 ** thewk(i), q(j,j,i), 1,_
                         qwk(j,j), 1)
    }
    #   main calculation
    call  dcopy (nobs, y, 1, ywk, 1)
    call  dcore (vmu, qwk, ldqr, nobs, n0, tol, ywk, job, limnla, nlawk,_
                 scrwk, varht, info, twk, work1)
    if (info != 0 )  return

    #   half the increment if no gain
    if ( scrold < scrwk ) {
        #   algorithm halts
        tmp = dabs (gwk1(idamax (nq, gwk1, 1)))
        if ( alph * tmp > - prec ) {
            info = -5
            return
        }
        alph = alph / 2.d0
        for (i=1;i<=nq;i=i+1)  thewk(i) = theta(i) + alph * gwk1(i)
        next
    }
    #   count for one iteration
    maxitwk = maxitwk - 1

    #   compute the gradient and the Hessian
    call  dcopy (n-2, qwk(n0+2,n0+1), ldqr+1, traux, 1)
    call  dcopy (n, qwk(n0+1,n0+1), ldqr+1, twk(2,1), 2)
    call  dcopy (n-1, qwk(n0+1,n0+2), ldqr+1, twk(1,2), 2)
    call  ddeev (vmu, nobs,_
                 q(n0+1,n0+1,1), ldqr, ldqc, n, nq, qwk(n0+2,n0+1),_
                 ldqr, traux, twk, ywk(n0+1),_                 
                 thewk, nlawk, scrwk, varht,_            # inputs
                 hes, nq, gra,_                          # outputs
                 hwk1, hwk2, gwk1, gwk2,_
                 kwk, n, work1, work2, c,_
                 info)

    #   get the active subset
    iwk = 0
    for (i=1;i<=nq;i=i+1) {
        if ( thewk(i) <= -25.d0 )  next
        iwk = iwk + 1
        call  dcopy (nq, hes(1,i), 1, hes(1,iwk), 1)
    }
    iwk = 0
    for (i=1;i<=nq;i=i+1) {
        if ( thewk(i) <= -25.d0 )  next
        iwk = iwk + 1
        call  dcopy (nq, hes(i,1), nq, hes(iwk,1), nq)
        gwk1(iwk) = gra(i)
        work2(iwk) = gra(i)
    }

    #   compute the Newton direction
    for (i=1;i<iwk;i=i+1)  
        call  dcopy (iwk-i, hes(i+1,i), 1, hes(i,i+1), nq)
    call  dmcdc (hes, nq, iwk, gwk2, pvtwk, info)
    call  dprmut (gwk1, iwk, pvtwk, 0)
    call  dposl (hes, nq, iwk, gwk1)
    call  dprmut (gwk1, iwk, pvtwk, 1)

    #   specify the stepsize
    alph = -1.d0

    #   set the update direction in the original index
    j = iwk
    for (i=nq;i>=1;i=i-1) {
        if ( thewk(i) <= -25.0 )  gwk1(i) = 0.d0
        else {
            gwk1(i) = gwk1(iwk)
            iwk = iwk - 1     
        }
    }
    call  dscal (nq, 1.d0/dlog(1.d1), gwk1, 1)
    tmp = dabs (gwk1(idamax (nq, gwk1, 1)))
    if ( tmp > 1.d0 )  call  dscal (nq, 1.d0/tmp, gwk1, 1)

    #   set thewk such that  nlawk = 0.d0
    for (i=1;i<=nq;i=i+1) {
        if ( thewk(i) <= -25.d0 )  next
        thewk(i) = thewk(i) - nlawk
    }
    call  dcopy (nq, thewk, 1, theta, 1)

    #   check convergence
    tmp = gra(idamax (nq, gra, 1)) ** 2
    if ( tmp < prec ** 2_                          #  zero gradient
        | scrold - scrwk < prec * (scrwk + 1.d0)_  #  small change
        & tmp < prec * (scrwk + 1.d0) ** 2 ) {     #  small gradient
        break
    }

    #   fail to converge
    if ( maxitwk < 1 ) {
        info = -4
        return
    }

    #   update records
    scrold = scrwk

    #   increment thewk
    for (i=1;i<=nq;i=i+1)  thewk(i) = thewk(i) + alph * gwk1(i)

    job = -1
    limnla(1) = -1.d0
    limnla(2) = 1.d0
}   #   the end of the main loop

#   compute the model to be returned
for (j=1;j<=nobs;j=j+1)  call  dset (nobs-j+1, 0.d0, qwk(j,j), 1)
for (i=1;i<=nq;i=i+1) {
    if ( theta(i) <= -25.d0 )  next
    for (j=1;j<=nobs;j=j+1)
        call  daxpy (nobs-j+1, 10.d0 ** theta(i), q(j,j,i), 1,_
                     qwk(j,j), 1)
}

call  dcopy (nobs, y, 1, ywk, 1)
call  dcore (vmu, qwk, ldqr, nobs, n0, tol, ywk, job, limnla, nlaht,_
             score, varht, info, twk, work1)
if (info != 0 )  return
call  dcoef (s, lds, nobs, n0, qraux, jpvt, ywk, qwk, ldqr, nlaht,_
             c, d, info, twk)

return
end

#....................................................................................
