
#::::::::::::::::::
#   cdennewton10
#::::::::::::::::::

subroutine  cdennewton10 (cd, nxis, q, nxi, rs, nobs, cntsum, cnt, intrs,
                          prec, maxiter, mchpr, jpvt, wk, info)

integer  nxis, nxi, nobs, cntsum, cnt(*), maxiter, jpvt(*), info
double precision  cd(*), q(nxi,*), rs(nobs,*), intrs(*), prec, mchpr, wk(*)

integer  iwt, imu, iv, icdnew, iwtnew, iwk

iwt = 1
imu = iwt + nobs
iv = imu + nxis
icdnew = iv + nxis*nxis
iwtnew = icdnew + nxis
iwk = iwtnew + nobs

call  cdennewton101 (cd, nxis, q, nxi, rs, nobs, cntsum, cnt, intrs,
                     prec, maxiter, mchpr, wk(iwt), wk(imu), wk(iv),
                     jpvt, wk(icdnew), wk(iwtnew), wk(iwk), info)

return
end


#:::::::::::::::::::
#   cdennewton101
#:::::::::::::::::::

subroutine  cdennewton101 (cd, nxis, q, nxi, rs, nobs, cntsum, cnt, intrs,
                           prec, maxiter, mchpr, wt, mu, v, jpvt, cdnew,
                           wtnew, wk, info)

integer  nxis, nxi, nobs, cntsum, cnt(*), maxiter, jpvt(*), info
double precision  cd(*), q(nxi,*), rs(nobs,*), intrs(*), prec, mchpr, wt(*),
                  mu(*), v(nxis,*), cdnew(*), wtnew(*), wk(*)

integer  i, j, k, iter, flag, rkv, idamax, infowk
double precision  tmp, ddot, dasum, wtsum, lkhd, mumax, wtsumnew, lkhdnew, disc, disc0

#   Initialization
info = 0
for (i=1;i<=nobs;i=i+1) {
    tmp = ddot (nxis, rs(i,1), nobs, cd, 1)
    wt(i) = dexp (-tmp)
    if (cntsum!=0)  wt(i) = wt(i) * dfloat (cnt(i))
}
wtsum = dasum (nobs, wt, 1)
lkhd = dlog (wtsum) + ddot (nxis, intrs, 1, cd, 1)
call  dsymv ('u', nxi, 1.d0, q, nxi, cd, 1, 0.d0, wk, 1)
lkhd = lkhd + ddot (nxi, cd, 1, wk, 1) / 2.d0
call  dscal (nobs, 1.d0/wtsum, wt, 1)
iter = 0
flag = 0
#   Iteration
repeat {
    iter = iter + 1
    # Calculate hessian and gradient
    for (i=1;i<=nxis;i=i+1)  mu(i) = ddot (nobs, wt, 1, rs(1,i), 1)
    for (i=1;i<=nxis;i=i+1) {
        for (j=i;j<=nxis;j=j+1) {
            v(i,j) = 0.d0
            for (k=1;k<=nobs;k=k+1)
                v(i,j) = v(i,j) + wt(k) * rs(k,i) * rs(k,j)
            v(i,j) = v(i,j) - mu(i) * mu(j)
            if (j<=nxi)  v(i,j) = v(i,j) + q(i,j)
        }
    }
    call  daxpy (nxis, -1.d0, intrs, 1, mu, 1)
    call  dsymv ('u', nxi, -1.d0, q, nxi, cd, 1, 1.d0, mu, 1)
    mumax = dabs(mu(idamax(nxis, mu, 1)))
    #   Cholesky factorization
    for (i=1;i<=nxis;i=i+1)  jpvt(i) = 0
    call  dchdc (v, nxis, nxis, wk, jpvt, 1, rkv)
    while (v(rkv,rkv)<v(1,1)*dsqrt(mchpr))  rkv = rkv - 1
    for (i=rkv+1;i<=nxis;i=i+1) {
        v(i,i) = v(1,1)
        call  dset (i-rkv-1, 0.d0, v(rkv+1,i), 1)
    }
    #   Update coefficients
    repeat {
        call  dcopy (nxis, mu, 1, cdnew, 1)
        call  dprmut (cdnew, nxis, jpvt, 0)
        call  dtrsl (v, nxis, nxis, cdnew, 11, infowk)
        call  dset (nxis-rkv, 0.d0, cdnew(rkv+1), 1)
        call  dtrsl (v, nxis, nxis, cdnew, 01, infowk)
        call  dprmut (cdnew, nxis, jpvt, 1)
        call  daxpy (nxis, 1.d0, cd, 1, cdnew, 1)
        for (i=1;i<=nobs;i=i+1) {
            tmp = ddot (nxis, rs(i,1), nobs, cdnew, 1)
            if (-tmp>3.d2) {
                flag = flag + 1
                break
            }
            wtnew(i) = dexp (-tmp)
            if (cntsum!=0)  wtnew(i) = wtnew(i) * dfloat (cnt(i))
        }
        wtsumnew = dasum (nobs, wtnew, 1)
        lkhdnew = dlog (wtsumnew) + ddot (nxis, intrs, 1, cdnew, 1)
        call  dsymv ('u', nxi, 1.d0, q, nxi, cdnew, 1, 0.d0, wk, 1)
        lkhdnew = lkhdnew + ddot (nxi, cdnew, 1, wk, 1) / 2.d0
        call  dscal (nobs, 1.d0/wtsumnew, wtnew, 1)
        #   Reset iteration with uniform starting value
        if (flag==1) {
            call  dset (nxis, 0.d0, cd, 1)
            if (cntsum!=0) {
                for (i=1;i<=nobs;i=i+1)  wt(i) = dfloat (cnt(i))
            }
            else  call  dset (nobs, 1.d0, wt, 1)
            wtsum = dasum (nobs, wt, 1)
            lkhd = dlog (wtsum)
            call  dscal (nobs, 1.d0/wtsum, wt, 1)
            iter = 0
            break
        }
        if (flag==3)  break
        if (lkhdnew-lkhd<1.d1*(1.d0+dabs(lkhd))*mchpr)  break
        call  dscal (nxis, .5d0, mu, 1)
        if (dabs(mu(idamax(nxis, mu, 1))/(1.d0+mumax))<1.d1*mchpr)  break
    }
    if (flag==1) {
        flag = 2
        next
    }
    if (flag==3) {
        info = 1
        return
    }
    #   Calculate convergence criterion
    disc = 0.d0
    for (i=1;i<=nobs;i=i+1)
        disc = dmax1 (disc, dabs(wt(i)-wtnew(i))/(1.d0+dabs(wt(i))))
    disc = dmax1 (disc, (mumax/(1.d0+dabs(lkhd)))**2)
    disc0 = dmax1 ((mumax/(1.d0+dabs(lkhd)))**2, dabs(lkhd-lkhdnew)/(1.d0+dabs(lkhd)))
    #   Set to new values
    call  dcopy (nxis, cdnew, 1, cd, 1)
    call  dcopy (nobs, wtnew, 1, wt, 1)
    lkhd = lkhdnew
    #   Check convergence
    if (disc0<prec)  break
    if (disc<prec)  break
    if (iter<maxiter)  next
    if (flag==0) {
        #   Reset iteration with uniform starting value
        call  dset (nxis, 0.d0, cd, 1)
        if (cntsum!=0) {
            for (i=1;i<=nobs;i=i+1)  wt(i) = dfloat (cnt(i))
        }
        else  call  dset (nobs, 1.d0, wt, 1)
        wtsum = dasum (nobs, wt, 1)
        lkhd = dlog (wtsum)
        call  dscal (nobs, 1.d0/wtsum, wt, 1)
        iter = 0
        flag = 2
    }
    else {
        info = 2
        break
    }
}
#   Calculate uncorrected v
for (i=1;i<=nxis;i=i+1) {
    for (j=i;j<=nxis;j=j+1) {
        v(i,j) = 0.d0
        for (k=1;k<=nobs;k=k+1)
            v(i,j) = v(i,j) + wt(k) * rs(k,i) * rs(k,j)
        if (j<=nxi)  v(i,j) = v(i,j) + q(i,j)
    }
}
#   Cholesky factorization
for (i=1;i<=nxis;i=i+1)  jpvt(i) = 0
call  dchdc (v, nxis, nxis, wk, jpvt, 1, rkv)
while (v(rkv,rkv)<v(1,1)*dsqrt(mchpr))  rkv = rkv - 1
for (i=rkv+1;i<=nxis;i=i+1) {
    v(i,i) = v(1,1)
    call  dset (i-rkv-1, 0.d0, v(rkv+1,i), 1)
}
# Calculate a
for (i=1;i<=nobs;i=i+1) {
    call  dcopy (nxis, rs(i,1), nobs, wk, 1)
    call  dprmut (wk, nxis, jpvt, 0)
    call  dtrsl (v, nxis, nxis, wk, 11, infowk)
    call  dset (nxis-rkv, 0.d0, wk(rkv+1), 1)
    wtnew(i) = wt(i) * ddot (nxis, wk, 1, wk, 1)
    if (cntsum!=0)  wtnew(i) = wtnew(i) / dfloat (cnt(i))
}
call  dcopy (nobs, wtnew, 1, wt, 1)

return
end
