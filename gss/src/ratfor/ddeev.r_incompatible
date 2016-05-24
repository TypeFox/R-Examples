
#:::::::::::
#   ddeev
#:::::::::::

subroutine  ddeev (vmu, nobs,_
                   q, ldqr, ldqc, n, nq, u, ldu, uaux, t, x,_     # inputs
                   theta, nlaht, score, varht,_              
                   hes, ldh, gra,_                                # outputs
                   hwk1, hwk2, gwk1, gwk2,_                       # work arrays
                   kwk, ldk, work1, work2, work3,_
                   info)

#  Acronym:  Double precision DErivative EValuation.

#  Purpose:  This routine calculates the gradient and the Hessian of
#      V(theta|lambda) or M(theta|lambda).
 
character*1       vmu
integer           nobs, ldqr, ldqc, n, nq, ldu, ldh, ldk, info
double precision  q(ldqr,ldqc,*), u(ldu,*), uaux(*), t(2,*), x(*),_
                  theta(*), nlaht, score, varht,_
                  hes(ldh,*), gra(*), hwk1(nq,*), hwk2(nq,*), gwk1(*), gwk2(*),_
                  kwk(ldk,ldk,*), work1(*), work2(*), work3(*)

#  On entry:
#      vmu        'v':  GCV criterion.
#                 'm':  GML criterion.
#                 'u':  unbiased risk estimate.
#      nobs       the number of observations.
#      q          F_{2}^{T} Q_{i} F_{2}, of size (n,n,nq).
#      n          the size of q.
#      nq         the number of Q_{i}'s.
#      u,uaux     Householder vectors of U, of size (n-1,n-2),
#                 where U^{T}DU is tridiagonal.
#      t          U^{T} (D-n\lambda I) U in packed form, of size (2,n).
#      x          U^{T}z = U^{T}F_{2}^{T}y, of size (n).
#      theta      the current log(theta) for the D matrix, of dimension (nq).
#      nlaht      the estimated  log10(n*lambda)  in the current model.
#      score      the minimum GCV/GML score found at (theta, nlaht).
#      varht      the variance estimate at (theta, nlaht).
 
#  On exit:
#      hes        Hessian at point (theta, nlaht), of size (nq,nq).
#      gra        gradient at point (theta, nlaht), of size (nq).
#      info        0 :  normal termination.
#                 -1 :  dimension error.
#                 -2 :  D !>= 0.
#                 -3 :  tuning parameters are out of scope.

#  Work arrays:
#      hwk1,2     of sizes at least (nq,nq).
#      gwk1,2     of sizes at least (nq).
#      kwk        of size at least (n,n,nq).
#      work1-3    of sizes at least (n).

#  Routines called directly:
#      Fortran -- dfloat
#      Blas    -- daxpy, dcopy, ddot, dscal
#      Blas2   -- dgemv
#      Linpack -- dpbfa, dpbsl, dqrsl
#      Rkpack  -- dqrslm
#      Other   -- dset

#  Written:  Chong Gu, Statistics, Purdue, latest version 12/29/91.

double precision  trc, det, dum, ddot
integer           i, j, m

info = 0
call  dset (nq, 0.d0, gra, 1)
call  dset (nq*nq, 0.d0, hes, 1)

#   check tuning parameters
if ( vmu != 'v' & vmu != 'm' & vmu != 'u' ) {
    info = -3
    return
}

#   check dimension
if ( nobs < n | ldqr < n | ldqc < n | nq <= 0 | ldu < n-1 | ldh < nq | ldk < n ) {
    info = -1
    return
}

#   compute  K_{i} = U^{T}(\theta_{i}Q_{i})U
for (i=2;i<=nq;i=i+1) {
#   from i=2 to nq
    if ( theta(i) <= -25.d0 ) next
    for (j=1;j<=n;j=j+1) {
        call  dcopy (n-j+1, q(j,j,i), 1, kwk(j,j,i), 1)
        call  dscal (n-j+1, 10.d0 ** theta(i), kwk(j,j,i), 1)
    }
    call  dqrslm (u, ldu, n-1, n-2, uaux, kwk(2,2,i), n, 0, info, work1)
    call  dqrsl (u, ldu, n-1, n-2, uaux, kwk(2,1,i), dum, kwk(2,1,i),_
                 dum, dum, dum, 01000, info)
}
#   compute K_{1} through the identity:  U^{T}(\sum K_{i})U = T
call  dcopy (n, t(2,1), 2, kwk(1,1,1), n+1)
call  dcopy (n-1, t(1,2), 2, kwk(2,1,1), n+1)
for (j=1;j<n-1;j=j+1)  call  dset (n-j-1, 0.d0, kwk(j+2,j,1), 1)
for (i=2;i<=nq;i=i+1) {
    if ( theta(i) <= -25.d0 )  next
    for (j=1;j<=n;j=j+1)
        call  daxpy (n-j+1, -1.d0, kwk(j,j,i), 1, kwk(j,j,1), 1)
}
#   fill the upper triangles of K_{i}
for (i=1;i<=nq;i=i+1) {
    if ( theta(i) <= -25.d0 )  next
    for (j=1;j<n;j=j+1)  call  dcopy (n-j, kwk(j+1,j,i), 1, kwk(j,j+1,i), n)
}

#   decompose the tridiagonal matrix  U^{T}DU
call  dset (n, 10.d0 ** nlaht, work1, 1)
call  daxpy (n, 1.d0, work1, 1, t(2,1), 2)
call  dpbfa (t, 2, n, 1, info)
             if ( info != 0 ) {
                 info = -2
                 return
             }

#   compute  T^{-1}K_{i}
for (i=1;i<=nq;i=i+1) {
    if ( theta(i) <= -25.d0 )  next
    for (j=1;j<=n;j=j+1)  call  dpbsl (t, 2, n, 1, kwk(1,j,i))
}

#::::::::::  Compute the gradient and the Hessian  ::::::::::

#   compute  -m x^{-T}T^{-m}K_{i}T^{-1}x:  m = 2('v') or 1('m')
call  dcopy (n, x, 1, work1, 1)
call  dpbsl (t, 2, n, 1, work1)
if ( vmu != 'm' ) {
    call  dcopy (n, work1, 1, work2, 1)
    call  dscal (n, 2.d0, work2, 1)
}
else  call  dcopy (n, x, 1, work2, 1)
for (i=1;i<=nq;i=i+1) {
    if ( theta(i) <= -25.d0 )  next
    call  dgemv ('t', n, n, 1.d0, kwk(1,1,i), n, work2, 1, 0.d0, work3, 1)
    gwk1(i) = - ddot (n, work1, 1, work3, 1)
}

#   compute  - tr[T^{-m}K_{i}]:  m = 2('v') or 1('m')
for (i=1;i<=nq;i=i+1) {
    gwk2(i) = 0.d0
    if ( theta(i) <= -25.d0 )  next
    for (j=1;j<=n;j=j+1) {
        if ( vmu != 'm' ) {
            call  dcopy (n, kwk(1,j,i), 1, work1, 1)
            call  dpbsl (t, 2, n, 1, work1)
            gwk2(i) = gwk2(i) - work1(j)
        }
        else  gwk2(i) = gwk2(i) - kwk(j,j,i)
    }
}
    
if ( vmu != 'm' ) {
    #   compute  2 x^{T}T^{-1} [K_{i}T^{-2}K_{j}+T^{-1}K_{i}T^{-1}K_{j}
    #                           +K_{i}T^{-1}K_{j}T^{-1}]T^{-1}x  for 'v'
    call  dcopy (n, x, 1, work1, 1)
    call  dpbsl (t, 2, n, 1, work1)
    for (i=1;i<=nq;i=i+1) {
        if ( theta(i) <= -25.d0 )  next
        call  dgemv ('n', n, n, 1.d0, kwk(1,1,i), n, work1, 1, 0.d0, work2, 1)
        for (j=1;j<=i;j=j+1) {
            if ( theta(j) <= -25.d0 )  next
            call  dgemv ('n', n, n, 1.d0, kwk(1,1,j), n, work1, 1, 0.d0, work3, 1)
            hwk1(i,j) = 2.d0 * ddot (n, work2, 1, work3, 1)
            call  dgemv ('t', n, n, 1.d0, kwk(1,1,j), n, work1, 1, 0.d0, work3, 1)
            hwk1(i,j) = hwk1(i,j) + 2.d0 * ddot (n, work2, 1, work3, 1)
        }
        call  dgemv ('t', n, n, 1.d0, kwk(1,1,i), n, work1, 1, 0.d0, work2, 1)
        for (j=1;j<=i;j=j+1) {
            if ( theta(j) <= -25.d0 )  next
            call  dgemv ('n', n, n, 1.d0, kwk(1,1,j), n, work1, 1, 0.d0, work3, 1)
            hwk1(i,j) = hwk1(i,j) + 2.d0 * ddot (n, work2, 1, work3, 1)
        }
    }
}
else {
    #   compute  2 x^{T} T^{-1}K_{i}T^{-1}K_{j}T^{-1}x  for 'm'
    call  dcopy (n, x, 1, work1, 1)
    call  dpbsl (t, 2, n, 1, work1)
    for (i=1;i<=nq;i=i+1) {
        if ( theta(i) <= -25.d0 )  next
        call  dgemv ('n', n, n, 1.d0, kwk(1,1,i), n, work1, 1, 0.d0, work2, 1)
        for (j=1;j<=i;j=j+1) {
            if ( theta(j) <= -25.d0 )  next
            call  dgemv ('t', n, n, 1.d0, kwk(1,1,j), n, x, 1, 0.d0, work3, 1)
            hwk1(i,j) = 2.d0 * ddot (n, work2, 1, work3, 1)
        }
    }
}
#   adjust diagonal
for (i=1;i<=nq;i=i+1) {
    if ( theta(i) <= -25.d0 )  next
    hwk1(i,i) = hwk1(i,i) + gwk1(i)
}
   
#   compute  m tr[T^{-m}K_{i}T^{-1}K_{j}]:  m = 2('v') or 1('m')
for (i=1;i<=nq;i=i+1) {
    if ( theta(i) <= -25.d0 )  next
    for (m=1;m<=i;m=m+1) {
        hwk2(i,m) = 0.d0
        if ( theta(m) <= -25.d0 )  next
        for (j=1;j<=n;j=j+1) {
            if ( vmu != 'm' ) {
                call  dcopy (n, kwk(1,j,m), 1, work1, 1)
                call  dpbsl (t, 2, n, 1, work1)
                hwk2(i,m) = hwk2(i,m) + 2.d0 * ddot (n, kwk(j,1,i), n, work1, 1)
            }
            else  hwk2(i,m) = hwk2(i,m) + ddot (n, kwk(j,1,i), n, kwk(1,j,m), 1)
        }
    }
}
#   adjust diagonal
for (i=1;i<=nq;i=i+1) {
    if ( theta(i) <= -25.d0 )  next
    hwk2(i,i) = hwk2(i,i) + gwk2(i)
}
    
#   compute the gradient
if ( vmu == 'v' ) {
    trc = dfloat (nobs) * 10.d0 ** (-nlaht) * varht / score
    for (i=1;i<=nq;i=i+1) {
        if ( theta(i) <= -25.d0 )  next
        gra(i) = gwk1(i) / trc / trc - 2.d0 * score * gwk2(i) / trc / dfloat(nobs)
    }
    call  dscal (nq, dfloat (nobs), gra, 1)
}
if ( vmu == 'u' ) {
    dum = 10.d0 ** nlaht
    for (i=1;i<=nq;i=i+1) {
        if ( theta(i) <= -25.d0 )  next
        gra(i) = dum * dum * gwk1(i) - 2.d0 * varht * dum * gwk2(i)
    }
    call  dscal (nq, 1.d0/dfloat (n), gra, 1)
}
if ( vmu == 'm' ) {
    det = 10.d0 ** (-nlaht) * varht / score
    for (i=1;i<=nq;i=i+1) {
        if ( theta(i) <= -25.d0 )  next
        gra(i) = gwk1(i) / det - dfloat (nobs) / dfloat (n) * score * gwk2(i)
    }
    call  dscal (nq, 1.d0 / dfloat (nobs), gra, 1)
}

#   compute the Hessian
if ( vmu == 'v' ) {
    for (i=1;i<=nq;i=i+1) {
        if ( theta(i) <= -25.d0 )  next
        for (j=1;j<=i;j=j+1) {
            if ( theta(j) <= -25.d0 )  next
            hes(i,j) = hwk1(i,j) / trc / trc - 2.d0 * gwk1(i) * gwk2(j) / trc ** 3_
                      - 2.d0 * gwk1(j) * gwk2(i) / trc ** 3 - 2.d0 * score * hwk2(i,j)_
                      / trc / dfloat (nobs) + 6.d0 * score * gwk2(i) * gwk2(j)_
                      / trc / trc / dfloat (nobs)
        }
        call  dscal (i, dfloat (nobs), hes(i,1), ldh)
    }
}
if ( vmu == 'u' ) {
    for (i=1;i<=nq;i=i+1) {
        if ( theta(i) <= -25.d0 )  next
        for (j=1;j<=i;j=j+1) {
            if ( theta(j) <= -25.d0 )  next
            hes(i,j) = dum * dum * hwk1(i,j) - 2.d0 * varht * dum * hwk2(i,j)
        }
        call  dscal (i, 1.d0/dfloat (n), hes(i,1), ldh)
    }
}
if ( vmu == 'm' ) {
    for (i=1;i<=nq;i=i+1) {
        if ( theta(i) <= -25.d0 )  next
        for (j=1;j<=i;j=j+1) {
            if ( theta(j) <= -25.d0 )  next
            hes(i,j) = hwk1(i,j) / det - gwk1(i) * gwk2(j) / det / dfloat (n)_
                      - gwk1(j) * gwk2(i) / det / dfloat (n) - dfloat (nobs)_
                      / dfloat (n) * score * hwk2(i,j) + dfloat (nobs)_
                      / dfloat (n) ** 2 * score * gwk2(i) * gwk2(j)
        }
        call  dscal (i, 1.d0 / dfloat (nobs), hes(i,1), ldh)
    }
}

return
end

#....................................................................................
