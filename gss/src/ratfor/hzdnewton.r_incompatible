
#:::::::::::::::
#   hzdnewton
#:::::::::::::::

subroutine  hzdnewton (cd, nxis, q, nxi, rs, nt, nobs, cntsum, cnt, qdrs, nqd,
                       qdwt, nx, prec, maxiter, mchpr, jpvt, wk, info)

integer  nxis, nxi, nt, nobs, cntsum, cnt(*), nqd, nx, maxiter, jpvt(*), info
double precision  cd(*), q(nxi,*), rs(nxis,*), qdrs(nqd,nxis,*), qdwt(nqd,*), prec, mchpr, wk(*)

integer  imrs, iwt, ifit, imu, imuwk, iv, ivwk, icdnew, iwtnew, ifitnew, iwk

imrs = 1
iwt = imrs + max0 (nxis, 2)
ifit = iwt + nqd*nx
imu = ifit + nt
imuwk = imu + nxis
iv = imuwk + nxis
ivwk = iv + nxis*nxis
icdnew = ivwk + nxis*nxis
iwtnew = icdnew + nxis
ifitnew = iwtnew + nqd*nx
iwk = ifitnew + nt

call  hzdnewton1 (cd, nxis, q, nxi, rs, nt, nobs, cntsum, cnt, qdrs, nqd, qdwt, nx,
                  prec, maxiter, mchpr, wk(imrs), wk(iwt), wk(ifit), wk(imu), wk(imuwk),
                  wk(iv), wk(ivwk), jpvt, wk(icdnew), wk(iwtnew), wk(ifitnew),
                  wk(iwk), info)

return
end


#::::::::::::::::
#   hzdnewton1
#::::::::::::::::

subroutine  hzdnewton1 (cd, nxis, q, nxi, rs, nt, nobs, cntsum, cnt, qdrs, nqd,
                        qdwt, nx, prec, maxiter, mchpr,
                        mrs, wt, fit, mu, muwk, v, vwk, jpvt, cdnew,
                        wtnew, fitnew, wk, info)

integer  nxis, nxi, nt, nobs, cntsum, cnt(*), nqd, nx, maxiter, jpvt(*), info
double precision  cd(*), q(nxi,*), rs(nxis,*), qdrs(nqd,nxis,*), qdwt(nqd,*),
                  prec, mchpr, mrs(*), wt(nqd,*), fit(*), mu(*), muwk(*),
                  v(nxis,*), vwk(nxis,*), cdnew(*), wtnew(nqd,*), fitnew(*), wk(*)

integer  i, j, k, kk, iter, flag, rkv, idamax, infowk
double precision  tmp, ddot, fitmean, dasum, lkhd, mumax, lkhdnew, disc, disc0, trc

#   Calculate constants
info = 0
for (i=1;i<=nxis;i=i+1) {
    mrs(i) = 0.d0
    for (j=1;j<=nt;j=j+1) {
        if (cntsum==0)  mrs(i) = mrs(i) + rs(i,j)
        else  mrs(i) = mrs(i) + rs(i,j) * dfloat (cnt(j))
    }
    mrs(i) = mrs(i) / dfloat (nobs)
}
#   Initialization
for (kk=1;kk<=nx;kk=kk+1) {
    for (i=1;i<=nqd;i=i+1)
        wt(i,kk) = qdwt(i,kk) * dexp (ddot (nxis, qdrs(i,1,kk), nqd, cd, 1))
}
fitmean = 0.d0
for (i=1;i<=nt;i=i+1) {
    tmp = ddot (nxis, rs(1,i), 1, cd, 1)
    fit(i) = dexp (tmp)
    if (cntsum!=0)  tmp = tmp * dfloat (cnt(i))
    fitmean = fitmean + tmp
}
fitmean = fitmean / dfloat (nobs) - dasum (nqd*nx, wt, 1)
call  dsymv ('u', nxi, 1.d0, q, nxi, cd, 1, 0.d0, wk, 1)
lkhd = ddot (nxi, cd, 1, wk, 1) / 2.d0 - fitmean
iter = 0
flag = 0
#   Iteration
repeat {
    iter = iter + 1
    # Calculate hessian and gradient
    call  dset (nxis, 0.d0, mu, 1)
    call  dset (nxis*nxis, 0.d0, v, 1)
    for (kk=1;kk<=nx;kk=kk+1) {
        for (i=1;i<=nxis;i=i+1) {
            muwk(i) = - ddot (nqd, wt(1,kk), 1, qdrs(1,i,kk), 1)
            for (j=i;j<=nxis;j=j+1) {
                vwk(i,j) = 0.d0
                for (k=1;k<=nqd;k=k+1)
                    vwk(i,j) = vwk(i,j) + wt(k,kk) * qdrs(k,i,kk) * qdrs(k,j,kk)
            }
        }
        call  daxpy (nxis, 1.d0, muwk, 1, mu, 1)
        call  daxpy (nxis*nxis, 1.d0, vwk, 1, v, 1)
    }
    for (i=1;i<=nxi;i=i+1) {
        for (j=i;j<=nxi;j=j+1)  v(i,j) = v(i,j) + q(i,j)
    }
    call  daxpy (nxis, 1.d0, mrs, 1, mu, 1)
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
        for (kk=1;kk<=nx;kk=kk+1) {
            for (i=1;i<=nqd;i=i+1) {
                tmp = ddot (nxis, qdrs(i,1,kk), nqd, cdnew, 1)
                if (tmp>3.d2) {
                    flag = flag + 1
                    break
                }
                wtnew(i,kk) = qdwt(i,kk) * dexp (tmp)
            }
            if ((flag==1)|(flag==3)) break
        }
        if ((flag==0)|(flag==2)) {
            fitmean = 0.d0
            for (i=1;i<=nt;i=i+1) {
                tmp = ddot (nxis, rs(1,i), 1, cdnew, 1)
                if (tmp>3.d2) {
                    flag = flag + 1
                    break
                }
                fitnew(i) = dexp (tmp)
                if (cntsum!=0)  tmp = tmp * dfloat (cnt(i))
                fitmean = fitmean + tmp
            }
            fitmean = fitmean / dfloat (nobs) - dasum (nqd*nx, wtnew, 1)
            call  dsymv ('u', nxi, 1.d0, q, nxi, cdnew, 1, 0.d0, wk, 1)
            lkhdnew = ddot (nxi, cdnew, 1, wk, 1) / 2.d0 - fitmean
        }
        #   Reset iteration with uniform starting value
        if (flag==1) {
            call  dset (nxis, 0.d0, cd, 1)
            call  dcopy (nqd*nx, qdwt, 1, wt, 1)
            fitmean = - dasum (nqd*nx, wt, 1)
            lkhd = - fitmean
            iter = 0
            break
        }
        if (flag==3)  break
        if (lkhdnew-lkhd<1.d1*(1.d0+dabs(lkhd))*mchpr)  break
        call  dscal (nxis, .5d0, mu, 1)
        if (dabs(mu(idamax(nxis, mu, 1))/mumax)<1.d1*mchpr)  break
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
    for (kk=1;kk<=nx;kk=kk+1) {
        for (i=1;i<=nqd;i=i+1)
            disc = dmax1 (disc, dabs(wt(i,kk)-wtnew(i,kk))/(1.d0+dabs(wt(i,kk))))
    }
    for (i=1;i<=nt;i=i+1)
        disc = dmax1 (disc, dabs(fit(i)-fitnew(i))/(1.d0+dabs(fit(i))))
    disc = dmax1 (disc, (mumax/(1.d0+dabs(lkhd)))**2)
    disc0 = dmax1 ((mumax/(1.d0+lkhd))**2, dabs(lkhd-lkhdnew)/(1.d0+dabs(lkhd)))
    #   Set to new values
    call  dcopy (nxis, cdnew, 1, cd, 1)
    call  dcopy (nqd*nx, wtnew, 1, wt, 1)
    call  dcopy (nt, fitnew, 1, fit, 1)
    lkhd = lkhdnew
    #   Check convergence
    if (disc0<prec)  break
    if (disc<prec)  break
    if (iter<maxiter)  next
    #   Reset iteration with uniform starting value
    if (flag==0) {
        call  dset (nxis, 0.d0, cd, 1)
        call  dcopy (nqd*nx, qdwt, 1, wt, 1)
        fitmean = - dasum (nqd*nx, wt, 1)
        lkhd = - fitmean
        iter = 0
        flag = 2
    }
    else {
        info = 2
        break
    }
}
#   Calculate proxy loss
for (i=1;i<=nt;i=i+1) {
    call  dprmut (rs(1,i), nxis, jpvt, 0)
    if (cntsum!=0)  call  dscal (nxis, dsqrt(dfloat(cnt(i))), rs(1,i), 1)
    call  dtrsl (v, nxis, nxis, rs(1,i), 11, infowk)
}
call  dprmut (mrs, nxis, jpvt, 0)
call  dtrsl (v, nxis, nxis, mrs, 11, infowk)
trc = ddot (nxis*nt, rs, 1, rs, 1) - dfloat (nobs) * ddot (nxis, mrs, 1, mrs, 1)
trc = trc / dfloat(nobs) / (dfloat(nobs)-1.d0)
mrs(1) = fitmean
mrs(2) = trc
for (kk=1;kk<=nx;kk=kk+1) {
    for (i=1;i<=nqd;i=i+1)
        wt(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cd, 1))
}

return
end
