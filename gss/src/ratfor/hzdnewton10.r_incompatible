
#:::::::::::::::::
#   hzdnewton10
#:::::::::::::::::

subroutine  hzdnewton10 (cd, nxis, q, nxi, rs, nt, nobs, cntsum, cnt,
                         intrs, rho, prec, maxiter, mchpr, jpvt, wk, info)

integer  nxis, nxi, nt, nobs, cntsum, cnt(*), maxiter, jpvt(*), info
double precision  cd(*), q(nxi,*), rs(nt,*), intrs(*), rho(*), prec, mchpr, wk(*)

integer  iwt, imu, iv, icdnew, iwtnew, iwk

iwt = 1
imu = iwt + nt
iv = imu + nxis
icdnew = iv + nxis*nxis
iwtnew = icdnew + nxis
iwk = iwtnew + nt

call  hzdnewton101 (cd, nxis, q, nxi, rs, nt, nobs, cntsum, cnt, intrs, rho,
                    prec, maxiter, mchpr, wk(iwt), wk(imu), wk(iv), jpvt,
                    wk(icdnew), wk(iwtnew), wk(iwk), info)

return
end


#::::::::::::::::::
#   hzdnewton101
#::::::::::::::::::

subroutine  hzdnewton101 (cd, nxis, q, nxi, rs, nt, nobs, cntsum, cnt,
                          intrs, rho, prec, maxiter, mchpr,
                          wt, mu, v, jpvt, cdnew, wtnew, wk, info)

integer  nxis, nxi, nt, nobs, cntsum, cnt(*), maxiter, jpvt(*), info
double precision  cd(*), q(nxi,*), rs(nt,*), intrs(*), rho(*), prec, mchpr,
                  wt(*), mu(*), v(nxis,*), cdnew(*), wtnew(*), wk(*)

integer  i, j, k, iter, flag, rkv, idamax, infowk
double precision  tmp, ddot, dasum, lkhd, mumax, lkhdnew, disc, disc0

#   Initialization
info = 0
for (i=1;i<=nt;i=i+1) {
    tmp = ddot (nxis, rs(i,1), nt, cd, 1)
    wt(i) = dexp (-tmp) * rho(i)
    if (cntsum!=0)  wt(i) = wt(i) * dfloat (cnt(i))
}
call  dscal (nt, 1/dfloat(nobs), wt, 1)
lkhd = dasum(nt, wt, 1) + ddot (nxis, intrs, 1, cd, 1)
call  dsymv ('u', nxi, 1.d0, q, nxi, cd, 1, 0.d0, wk, 1)
lkhd = lkhd + ddot (nxi, cd, 1, wk, 1) / 2.d0
iter = 0
flag = 0
#   Iteration
repeat {
    iter = iter + 1
    # Calculate hessian and gradient
    for (i=1;i<=nxis;i=i+1) {
        mu(i) = ddot (nt, wt, 1, rs(1,i), 1)
        for (j=i;j<=nxis;j=j+1) {
            v(i,j) = 0.d0
            for (k=1;k<=nt;k=k+1)
                v(i,j) = v(i,j) + wt(k) * rs(k,i) * rs(k,j)
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
        for (i=1;i<=nt;i=i+1) {
            tmp = ddot (nxis, rs(i,1), nt, cdnew, 1)
            if (-tmp>3.d2) {
                flag = flag + 1
                break
            }
            wtnew(i) = dexp (-tmp) * rho(i)
            if (cntsum!=0)  wtnew(i) = wtnew(i) * dfloat (cnt(i))
        }
        call  dscal (nt, 1/dfloat(nobs), wtnew, 1)
        lkhdnew = dasum(nt, wtnew, 1) + ddot (nxis, intrs, 1, cdnew, 1)
        call  dsymv ('u', nxi, 1.d0, q, nxi, cdnew, 1, 0.d0, wk, 1)
        lkhdnew = lkhdnew + ddot (nxi, cdnew, 1, wk, 1) / 2.d0
        #   Reset iteration with uniform starting value
        if (flag==1) {
            call  dset (nxis, 0.d0, cd, 1)
            for (i=1;i<=nt;i=i+1) {
                wt(i) = rho(i)
                if (cntsum!=0)  wt(i) = wt(i) * dfloat (cnt(i))
            }
            call  dscal (nt, 1/dfloat(nobs), wt, 1)
            lkhd = dasum (nt, wt, 1)
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
    for (i=1;i<=nt;i=i+1)
        disc = dmax1 (disc, dabs(wt(i)-wtnew(i))/(1.d0+dabs(wt(i))))
    disc = dmax1 (disc, (mumax/(1.d0+dabs(lkhd)))**2)
    disc0 = dmax1 ((mumax/(1.d0+dabs(lkhd)))**2, dabs(lkhd-lkhdnew)/(1.d0+dabs(lkhd)))
    #   Set to new values
    call  dcopy (nxis, cdnew, 1, cd, 1)
    call  dcopy (nt, wtnew, 1, wt, 1)
    lkhd = lkhdnew
    #   Check convergence
    if (disc0<prec)  break
    if (disc<prec)  break
    if (iter<maxiter)  next
    if (flag==0) {
        #   Reset iteration with uniform starting value
        call  dset (nxis, 0.d0, cd, 1)
        for (i=1;i<=nt;i=i+1) {
            wt(i) = rho(i)
            if (cntsum!=0)  wt(i) = wt(i) * dfloat (cnt(i))
        }
        call  dscal (nt, 1/dfloat(nobs), wt, 1)
        lkhd = dasum (nt, wt, 1)
        iter = 0
        flag = 2
    }
    else {
        info = 2
        break
    }
}
#   Calculate proxy loss
lkhd = dasum (nt, wt, 1) + ddot (nxis, intrs, 1, cd, 1)
tmp = 0.d0
disc = 0.d0
for (i=1;i<=nt;i=i+1) {
    call  dcopy (nxis, rs(i,1), nt, wk, 1)
    call  dprmut (wk, nxis, jpvt, 0)
    call  dtrsl (v, nxis, nxis, wk, 11, infowk)
    call  dset (nxis-rkv, 0.d0, wk(rkv+1), 1)
    wtnew(i) = wt(i) * ddot (nxis, wk, 1, wk, 1)
    if (cntsum!=0)  wtnew(i) = wtnew(i) / dfloat (cnt(i))
    tmp = tmp + wt(i) * (dexp (wtnew(i)/(1.d0-wtnew(i))) - 1.d0)
#    tmp = tmp + wt(i) * (1.d0/(1.d0-wtnew(i))**2-1.d0)/2.d0
    if (cntsum!=0) disc = disc + dfloat(cnt(i)) * wtnew(i)/(1.d0-wtnew(i))
    else  disc = disc + wtnew(i)/(1.d0-wtnew(i))
}
wt(1) = lkhd
wt(2) = tmp
wt(3) = disc/dfloat(nobs)

return
end


#:::::::::::::::
#   hzdaux101
#:::::::::::::::

subroutine  hzdaux101 (cd, nxis, q, nxi, rs, nt, rho, mchpr, v, jpvt)

integer  nxis, nxi, nt, jpvt(*)
double precision  cd(*), q(nxi,*), rs(nt,*), rho(*), mchpr, v(nxis,*)

integer  i, j, k, rkv
double precision  tmp, ddot

#   Initialization
for (i=1;i<=nt;i=i+1) {
    tmp = ddot (nxis, rs(i,1), nt, cd, 1)
    rho(i) = dexp (-tmp) * rho(i)
}
#   H matrix
for (i=1;i<=nxis;i=i+1) {
    for (j=i;j<=nxis;j=j+1) {
        v(i,j) = 0.d0
        for (k=1;k<=nt;k=k+1)
            v(i,j) = v(i,j) + rho(k) * rs(k,i) * rs(k,j)
        if (j<=nxi)  v(i,j) = v(i,j) + q(i,j)
    }
}
#   Cholesky factorization
for (i=1;i<=nxis;i=i+1)  jpvt(i) = 0
call  dchdc (v, nxis, nxis, cd, jpvt, 1, rkv)
while (v(rkv,rkv)<v(1,1)*dsqrt(mchpr))  rkv = rkv - 1
for (i=rkv+1;i<=nxis;i=i+1) {
    v(i,i) = v(1,1)
    call  dset (i-rkv-1, 0.d0, v(rkv+1,i), 1)
}

return
end
