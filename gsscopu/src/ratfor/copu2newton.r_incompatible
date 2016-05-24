
#:::::::::::::::::
#   copu2newton
#:::::::::::::::::

subroutine  copu2newton (cd, nxis, q, nxi, rs0, n0, cntsum0, cnt0, qdrs, nqd,
                         qdrs1, wt1, n1, cntsum1, cnt1, qdrs2, wt2, n2, cntsum2,
                         cnt2, wt3, n3, cntsum3, cnt3, nt, twt, qdwt, tind,
                         prec, maxiter, mchpr, jpvt, wk, info)

integer  nxis, nxi, n0, cntsum0, cnt0(*), nqd, n1, cntsum1, cnt1(*), n2, cntsum2,
         cnt2(*), n3, cntsum3, cnt3(*), nt, tind(*), maxiter, jpvt(*), info
double precision  cd(*), q(nxi,*), rs0(nxis,*), qdrs(nxis,*), qdrs1(nqd,nxis,*),
                  wt1(nqd,*), qdrs2(nqd,nxis,*), wt2(nqd,*), wt3(nqd,2,*), twt(*),
                  qdwt(nqd,2,*), prec, mchpr, wk(*)

integer  imrs, imrs2, ieta, ieta1, ieta2, imu, imuwk, iv, ivwk, icdnew, imut, iwk

imrs = 1
imrs2 = imrs + nxis
ieta = imrs2 + nxis
ieta1 = ieta + nqd*nqd
ieta2 = ieta1 + nqd*n1
imu = ieta2 + nqd*n2
imuwk = imu + nxis
iv = imuwk + nxis
ivwk = iv + nxis*nxis
icdnew = ivwk + nxis*nxis
imut = icdnew + nxis
iwk = imut + nxis*nt

call  copu2newton1 (cd, nxis, q, nxi, rs0, n0, cntsum0, cnt0, qdrs, nqd, qdrs1,
                    wt1, n1, cntsum1, cnt1, qdrs2, wt2, n2, cntsum2, cnt2, wt3,
                    n3, cntsum3, cnt3, nt, twt, qdwt, tind, prec, maxiter, mchpr,
                    wk(imrs), wk(imrs2), wk(ieta), wk(ieta1), wk(ieta2), wk(imu),
                    wk(imuwk), wk(iv), wk(ivwk), jpvt, wk(icdnew), wk(imut),
                    wk(iwk), info)

return
end


#::::::::::::::::::
#   copu2newton1
#::::::::::::::::::

subroutine  copu2newton1 (cd, nxis, q, nxi, rs0, n0, cntsum0, cnt0, qdrs, nqd,
                          qdrs1, wt1, n1, cntsum1, cnt1, qdrs2, wt2, n2, cntsum2,
                          cnt2, wt3, n3, cntsum3, cnt3, nt, twt, qdwt, tind,
                          prec, maxiter, mchpr, mrs, mrs2, eta, eta1, eta2, mu,
                          muwk, v, vwk, jpvt, cdnew, mut, wk, info)

integer  nxis, nxi, n0, cntsum0, cnt0(*), nqd, n1, cntsum1, cnt1(*), n2, cntsum2,
         cnt2(*), n3, cntsum3, cnt3(*), nt, tind(*), maxiter, jpvt(*), info
double precision  cd(*), q(nxi,*), rs0(nxis,*), qdrs(nxis,*), qdrs1(nqd,nxis,*),
                  wt1(nqd,*), qdrs2(nqd,nxis,*), wt2(nqd,*), wt3(nqd,2,*), twt(*),
                  qdwt(nqd,2,*), prec, mchpr, mrs(*), mrs2(*), eta(*), eta1(nqd,*),
                  eta2(nqd,*), mu(*), muwk(*), v(nxis,*), vwk(nxis,*), cdnew(*),
                  mut(nxis,*), wk(*)

integer  nobs, i, j, k, kk, iter, flag, rkv, idamax, infowk
double precision  tmp, ddot, dasum, lkhd, mumax, lkhdnew, disc, disc0, trc

#   Calculate constants
info = 0
if (cntsum0==0)  nobs = n0
else  nobs = cntsum0
if (cntsum1==0)  nobs = nobs + n1
else  nobs = nobs + cntsum1
if (cntsum2==0)  nobs = nobs + n2
else  nobs = nobs + cntsum2
if (cntsum3==0)  nobs = nobs + n3
else  nobs = nobs + cntsum3
for (i=1;i<=nxis;i=i+1) {
    mrs(i) = 0.d0
    for (j=1;j<=n0;j=j+1) {
        if (cntsum0==0)  mrs(i) = mrs(i) + rs0(i,j)
        else  mrs(i) = mrs(i) + rs0(i,j) * dfloat (cnt0(j))
    }
}
#   Initialization
for (i=1;i<=nqd*nqd;i=i+1)  eta(i) = dexp (ddot (nxis, qdrs(1,i), 1, cd, 1))
lkhd = 0.d0
for (i=1;i<=n0;i=i+1) {
    tmp = ddot (nxis, rs0(1,i), 1, cd, 1)
    if (cntsum0!=0)  tmp = tmp * dfloat (cnt0(i))
    lkhd = lkhd - tmp
}
for (i=1;i<=n1;i=i+1) {
    for (j=1;j<=nqd;j=j+1) {
        eta1(j,i) = dexp (ddot (nxis, qdrs1(j,1,i), nqd, cd, 1)) * wt1(j,i)
    }
    if (cntsum1==0)  lkhd = lkhd - dlog (dasum (nqd, eta1(1,i), 1))
    else  lkhd = lkhd - dlog (dasum (nqd, eta1(1,i), 1)) * dfloat (cnt1(i))
}
for (i=1;i<=n2;i=i+1) {
    for (j=1;j<=nqd;j=j+1) {
        eta2(j,i) = dexp (ddot (nxis, qdrs2(j,1,i), nqd, cd, 1)) * wt2(j,i)
    }
    if (cntsum2==0)  lkhd = lkhd - dlog (dasum (nqd, eta2(1,i), 1))
    else  lkhd = lkhd - dlog (dasum (nqd, eta2(1,i), 1)) * dfloat (cnt2(i))
}
for (i=1;i<=n3;i=i+1) {
    tmp = 0.d0
    for (j=1;j<=nqd;j=j+1) {
        tmp = tmp + ddot (nqd, eta((j-1)*nqd+1), 1, wt3(1,1,i), 1) * wt3(j,2,i)
    }
    if (cntsum3==0)  lkhd = lkhd - dlog (tmp)
    else  lkhd = lkhd - dlog (tmp) * dfloat (cnt3(i))
}
lkhd = lkhd / dfloat (nobs)
for (i=1;i<=nt;i=i+1) {
    tmp = 0.d0
    for (j=1;j<=nqd;j=j+1) {
        tmp = tmp + ddot (nqd, eta((j-1)*nqd+1), 1, qdwt(1,1,i), 1) * qdwt(j,2,i)
    }
    lkhd = lkhd + dlog (tmp) * twt(i)
}
call  dsymv ('u', nxi, 1.d0, q, nxi, cd, 1, 0.d0, muwk, 1)
lkhd = lkhd + ddot (nxi, cd, 1, muwk, 1) / 2.d0
iter = 0
flag = 0
#   Iteration
repeat {
    iter = iter + 1
    #   Calculate hessian and gradient
    call  dcopy (nxis, mrs, 1, mu, 1)
    call  dset(nxis*nxis, 0.d0, v, 1)
    for (i=1;i<=n1;i=i+1) {
        tmp = dasum (nqd, eta1(1,i), 1)
        for (j=1;j<=nxis;j=j+1)
            muwk(j) = ddot (nqd, eta1(1,i), 1, qdrs1(1,j,i), 1) / tmp
        for (j=1;j<=nxis;j=j+1) {
            for (k=j;k<=nxis;k=k+1) {
                vwk(j,k) = 0.d0
                for (kk=1;kk<=nqd;kk=kk+1)
                    vwk(j,k) = vwk(j,k) + eta1(kk,i)*qdrs1(kk,j,i)*qdrs1(kk,k,i)
                vwk(j,k) = vwk(j,k) / tmp - muwk(j) * muwk(k)
            }
        }
        if (cntsum1==0) {
            call  daxpy (nxis, 1.d0, muwk, 1, mu, 1)
            call  daxpy (nxis*nxis, -1.d0, vwk, 1, v, 1)
        }
        else {
            call  daxpy (nxis, dfloat (cnt1(i)), muwk, 1, mu, 1)
            call  daxpy (nxis*nxis, -dfloat (cnt1(i)), vwk, 1, v, 1)
        }
    }
    for (i=1;i<=n2;i=i+1) {
        tmp = dasum (nqd, eta2(1,i), 1)
        for (j=1;j<=nxis;j=j+1)
            muwk(j) = ddot (nqd, eta2(1,i), 1, qdrs2(1,j,i), 1) / tmp
        for (j=1;j<=nxis;j=j+1) {
            for (k=j;k<=nxis;k=k+1) {
                vwk(j,k) = 0.d0
                for (kk=1;kk<=nqd;kk=kk+1)
                    vwk(j,k) = vwk(j,k) + eta2(kk,i)*qdrs2(kk,j,i)*qdrs2(kk,k,i)
                vwk(j,k) = vwk(j,k) / tmp - muwk(j) * muwk(k)
            }
        }
        if (cntsum2==0) {
            call  daxpy (nxis, 1.d0, muwk, 1, mu, 1)
            call  daxpy (nxis*nxis, -1.d0, vwk, 1, v, 1)
        }
        else {
            call  daxpy (nxis, dfloat (cnt2(i)), muwk, 1, mu, 1)
            call  daxpy (nxis*nxis, -dfloat (cnt2(i)), vwk, 1, v, 1)
        }
    }
    for (i=1;i<=n3;i=i+1) {
        for (j=1;j<=nqd;j=j+1) {
            for (k=1;k<=nqd;k=k+1)
                wk((k-1)*nqd+j) = eta((k-1)*nqd+j) * wt3(j,1,i) * wt3(k,2,i)
        }
        tmp = dasum (nqd*nqd, wk, 1)
        for (j=1;j<=nxis;j=j+1)
            muwk(j) = ddot (nqd*nqd, wk, 1, qdrs(j,1), nxis) / tmp
        for (j=1;j<=nxis;j=j+1) {
            for (k=j;k<=nxis;k=k+1) {
                vwk(j,k) = 0.d0
                for (kk=1;kk<=nqd*nqd;kk=kk+1)
                    vwk(j,k) = vwk(j,k) + wk(kk)*qdrs(j,kk)*qdrs(k,kk)
                vwk(j,k) = vwk(j,k) / tmp - muwk(j) * muwk(k)
            }
        }
        if (cntsum3==0) {
            call  daxpy (nxis, 1.d0, muwk, 1, mu, 1)
            call  daxpy (nxis*nxis, -1.d0, vwk, 1, v, 1)
        }
        else {
            call  daxpy (nxis, dfloat (cnt3(i)), muwk, 1, mu, 1)
            call  daxpy (nxis*nxis, -dfloat (cnt3(i)), vwk, 1, v, 1)
        }
    }
    call  dscal (nxis, 1.d0/dfloat(nobs), mu, 1)
    call  dscal (nxis*nxis, 1.d0/dfloat(nobs), v, 1)
    for (i=1;i<=nt;i=i+1) {
        for (j=1;j<=nqd;j=j+1) {
            for (k=1;k<=nqd;k=k+1)
                wk((k-1)*nqd+j) = eta((k-1)*nqd+j) * qdwt(j,1,i) * qdwt(k,2,i)
        }
        tmp = dasum (nqd*nqd, wk, 1)
        for (j=1;j<=nxis;j=j+1)
            muwk(j) = ddot (nqd*nqd, wk, 1, qdrs(j,1), nxis) / tmp
        for (j=1;j<=nxis;j=j+1) {
            for (k=j;k<=nxis;k=k+1) {
                vwk(j,k) = 0.d0
                for (kk=1;kk<=nqd*nqd;kk=kk+1)
                    vwk(j,k) = vwk(j,k) + wk(kk)*qdrs(j,kk)*qdrs(k,kk)
                vwk(j,k) = vwk(j,k) / tmp - muwk(j) * muwk(k)
            }
        }
        call  dcopy (nxis, muwk, 1, mut(1,i), 1)
        call  daxpy (nxis, -twt(i), muwk, 1, mu, 1)
        call  daxpy (nxis*nxis, twt(i), vwk, 1, v, 1)
    }
    call  dcopy (nxis, mu, 1, mrs2, 1)
    call  dsymv ('u', nxi, -1.d0, q, nxi, cd, 1, 1.d0, mu, 1)
    for (i=1;i<=nxi;i=i+1) {
        for (j=i;j<=nxi;j=j+1)  v(i,j) = v(i,j) + q(i,j)
    }
    mumax = dabs(mu(idamax(nxis, mu, 1)))
    #   Cholesky factorization
    for (i=1;i<=nxis;i=i+1)  jpvt(i) = 0
#    call  dmcdc (v, nxis, nxis, muwk, jpvt, infowk)
    call  dchdc (v, nxis, nxis, muwk, jpvt, 1, rkv)
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
        for (i=1;i<=nqd*nqd;i=i+1) {
            tmp = ddot (nxis, qdrs(1,i), 1, cdnew, 1)
            if (tmp>3.d2) {
                flag = flag + 1
                break
            }
            wk(i) = dexp (tmp)
        }
        lkhdnew = 0.d0
        for (i=1;i<=n0;i=i+1) {
            tmp = ddot (nxis, rs0(1,i), 1, cdnew, 1)
            if (cntsum0!=0)  tmp = tmp * dfloat (cnt0(i))
            lkhdnew = lkhdnew - tmp
        }
        for (i=1;i<=n1;i=i+1) {
            for (j=1;j<=nqd;j=j+1) {
                eta1(j,i) = dexp (ddot (nxis, qdrs1(j,1,i), nqd, cdnew, 1)) * wt1(j,i)
            }
            if (cntsum1==0)  lkhdnew = lkhdnew - dlog (dasum (nqd, eta1(1,i), 1))
            else  lkhdnew = lkhdnew - dlog (dasum (nqd, eta1(1,i), 1)) * dfloat (cnt1(i))
        }
        for (i=1;i<=n2;i=i+1) {
            for (j=1;j<=nqd;j=j+1) {
                eta2(j,i) = dexp (ddot (nxis, qdrs2(j,1,i), nqd, cdnew, 1)) * wt2(j,i)
            }
            if (cntsum2==0)  lkhdnew = lkhdnew - dlog (dasum (nqd, eta2(1,i), 1))
            else  lkhdnew = lkhdnew - dlog (dasum (nqd, eta2(1,i), 1)) * dfloat (cnt2(i))
        }
        for (i=1;i<=n3;i=i+1) {
            tmp = 0.d0
            for (j=1;j<=nqd;j=j+1) {
                tmp = tmp + ddot (nqd, wk((j-1)*nqd+1), 1, wt3(1,1,i), 1) * wt3(j,2,i)
            }
            if (cntsum3==0)  lkhdnew = lkhdnew - dlog (tmp)
            else  lkhdnew = lkhdnew - dlog (tmp) * dfloat (cnt3(i))
        }
        lkhdnew = lkhdnew / dfloat (nobs)
        for (i=1;i<=nt;i=i+1) {
            tmp = 0.d0
            for (j=1;j<=nqd;j=j+1) {
                tmp = tmp + ddot (nqd, wk((j-1)*nqd+1), 1, qdwt(1,1,i), 1) * qdwt(j,2,i)
            }
            lkhdnew = lkhdnew + dlog (tmp) * twt(i)
        }
        call  dsymv ('u', nxi, 1.d0, q, nxi, cdnew, 1, 0.d0, muwk, 1)
        lkhdnew = lkhdnew + ddot (nxi, cdnew, 1, muwk, 1) / 2.d0
        #   Reset iteration with uniform starting value
        if (flag==1) {
            call  dset (nxis, 0.d0, cd, 1)
            call  dset (nqd*nqd, 1.d0, eta, 1)
            call  dcopy (nqd*n1, wt1, 1, eta1, 1)
            call  dcopy (nqd*n2, wt2, 1, eta2, 1)
            lkhd = 0.d0
            for (i=1;i<=n1;i=i+1) {
                if (cntsum1==0)  lkhd = lkhd - dlog (dasum (nqd, eta1(1,i), 1))
                else  lkhd = lkhd - dlog (dasum (nqd, eta1(1,i), 1)) * dfloat (cnt1(i))
            }
            for (i=1;i<=n2;i=i+1) {
                if (cntsum2==0)  lkhd = lkhd - dlog (dasum (nqd, eta2(1,i), 1))
                else  lkhd = lkhd - dlog (dasum (nqd, eta2(1,i), 1)) * dfloat (cnt2(i))
            }
            for (i=1;i<=n3;i=i+1) {
                tmp = dasum (nqd, wt3(1,1,i), 1) * dasum (nqd, wt3(1,2,i), 1)
                if (cntsum3==0)  lkhdnew = lkhdnew - dlog (tmp)
                else  lkhdnew = lkhdnew - dlog (tmp) * dfloat (cnt3(i))
            }
            lkhd = lkhd / dfloat (nobs)
            for (i=1;i<=nt;i=i+1) {
                tmp = dasum (nqd, qdwt(1,1,i), 1) * dasum (nqd, qdwt(1,2,i), 1)
                lkhd = lkhd + dlog (tmp) * twt(i)
            }
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
    for (i=1;i<=nqd*nqd;i=i+1) {
        disc = dmax1 (disc, dabs(eta(i)-wk(i))/(1.d0+dabs(eta(i))))
    }
    disc = dmax1 (disc, (mumax/(1.d0+dabs(lkhd)))**2)
    disc0 = dmax1 ((mumax/(1.d0+lkhd))**2, dabs(lkhd-lkhdnew)/(1.d0+dabs(lkhd)))
    #   Set to new values
    call  dcopy (nxis, cdnew, 1, cd, 1)
    call  dcopy (nqd*nqd, wk, 1, eta, 1)
    lkhd = lkhdnew
    #   Check convergence
    if (disc0<prec)  break
    if (disc<prec)  break
    if (iter<maxiter)  next
    if (flag==0) {
        #   Reset iteration with uniform starting value
        call  dset (nxis, 0.d0, cd, 1)
        call  dset (nqd*nqd, 1.d0, eta, 1)
        call  dcopy (nqd*n1, wt1, 1, eta1, 1)
        call  dcopy (nqd*n2, wt2, 1, eta2, 1)
        lkhd = 0.d0
        for (i=1;i<=n1;i=i+1) {
            if (cntsum1==0)  lkhd = lkhd - dlog (dasum (nqd, eta1(1,i), 1))
            else  lkhd = lkhd - dlog (dasum (nqd, eta1(1,i), 1)) * dfloat (cnt1(i))
        }
        for (i=1;i<=n2;i=i+1) {
            if (cntsum2==0)  lkhd = lkhd - dlog (dasum (nqd, eta2(1,i), 1))
            else  lkhd = lkhd - dlog (dasum (nqd, eta2(1,i), 1)) * dfloat (cnt2(i))
        }
        for (i=1;i<=n3;i=i+1) {
            tmp = dasum (nqd, wt3(1,1,i), 1) * dasum (nqd, wt3(1,2,i), 1)
            if (cntsum3==0)  lkhdnew = lkhdnew - dlog (tmp)
            else  lkhdnew = lkhdnew - dlog (tmp) * dfloat (cnt3(i))
        }
        lkhd = lkhd / dfloat (nobs)
        for (i=1;i<=nt;i=i+1) {
            tmp = dasum (nqd, qdwt(1,1,i), 1) * dasum (nqd, qdwt(1,2,i), 1)
            lkhd = lkhd + dlog (tmp) * twt(i)
        }
        iter = 0
        flag = 2
    }
    else {
        info = 2
        break
    }
}
#   Calculate proxy loss
trc = 0.d0
for (i=1;i<=n0;i=i+1) {
    call  dcopy (nxis, rs0(1,i), 1, muwk, 1)
    if (nt>1)  call  daxpy (nxis, -1.d0, mut(1,tind(i)), 1, muwk, 1)
    else  call  daxpy (nxis, -1.d0, mut, 1, muwk, 1)
    call  daxpy (nxis, -1.d0, mrs2, 1, muwk, 1)
    call  dprmut (muwk, nxis, jpvt, 0)
    if (cntsum0!=0)  call  dscal (nxis, dsqrt(dfloat(cnt0(i))), muwk, 1)
    call  dtrsl (v, nxis, nxis, muwk, 11, infowk)
    if (nxis-rkv>0)  call  dset (nxis-rkv, 0.d0, muwk(rkv+1), 1)
    trc = trc + ddot (nxis, muwk, 1, muwk, 1)
}
for (i=1;i<=n1;i=i+1) {
    tmp = dasum (nqd, eta1(1,i), 1)
    for (j=1;j<=nxis;j=j+1)
        muwk(j) = ddot (nqd, eta1(1,i), 1, qdrs1(1,j,i), 1) / tmp
    if (nt>1)  call  daxpy (nxis, -1.d0, mut(1,tind(n0+i)), 1, muwk, 1)
    else  call  daxpy (nxis, -1.d0, mut, 1, muwk, 1)
    call  daxpy (nxis, -1.d0, msr2, 1, muwk, 1)
    call  dprmut (muwk, nxis, jpvt, 0)
    if (cntsum1!=0)  call  dscal (nxis, dsqrt(dfloat(cnt1(i))), muwk, 1)
    call  dtrsl (v, nxis, nxis, muwk, 11, infowk)
    if (nxis-rkv>0)  call  dset (nxis-rkv, 0.d0, muwk(rkv+1), 1)
    trc = trc + ddot (nxis, muwk, 1, muwk, 1)
}
for (i=1;i<=n2;i=i+1) {
    tmp = dasum (nqd, eta2(1,i), 1)
    for (j=1;j<=nxis;j=j+1)
        muwk(j) = ddot (nqd, eta2(1,i), 1, qdrs2(1,j,i), 1) / tmp
    if (nt>1)  call  daxpy (nxis, -1.d0, mut(1,tind(n0+n1+i)), 1, muwk, 1)
    else  call  daxpy (nxis, -1.d0, mut, 1, muwk, 1)
    call  daxpy (nxis, -1.d0, mrs2, 1, muwk, 1)
    call  dprmut (muwk, nxis, jpvt, 0)
    if (cntsum2!=0)  call  dscal (nxis, dsqrt(dfloat(cnt2(i))), muwk, 1)
    call  dtrsl (v, nxis, nxis, muwk, 11, infowk)
    if (nxis-rkv>0)  call  dset (nxis-rkv, 0.d0, muwk(rkv+1), 1)
    trc = trc + ddot (nxis, muwk, 1, muwk, 1)
}
for (i=1;i<=n3;i=i+1) {
    for (j=1;j<=nqd;j=j+1) {
        for (k=1;k<=nqd;k=k+1)
            wk((k-1)*nqd+j) = eta((k-1)*nqd+j) * wt3(j,1,i) * wt3(k,2,i)
    }
    tmp = dasum (nqd*nqd, wk, 1)
    for (j=1;j<=nxis;j=j+1)
        muwk(j) = ddot (nqd*nqd, wk, 1, qdrs(j,1), nxis) / tmp
    if (nt>1)  call  daxpy (nxis, -1.d0, mut(1,tind(n0+n1+n2+i)), 1, muwk, 1)
    else  call  daxpy (nxis, -1.d0, mut, 1, muwk, 1)
    call  daxpy (nxis, -1.d0, mrs2, 1, muwk, 1)
    call  dprmut (muwk, nxis, jpvt, 0)
    if (cntsum3!=0)  call  dscal (nxis, dsqrt(dfloat(cnt3(i))), muwk, 1)
    call  dtrsl (v, nxis, nxis, muwk, 11, infowk)
    if (nxis-rkv>0)  call  dset (nxis-rkv, 0.d0, muwk(rkv+1), 1)
    trc = trc + ddot (nxis, muwk, 1, muwk, 1)
}
trc = trc / dfloat(nobs) / (dfloat(nobs)-1.d0)
call  dsymv ('u', nxi, 1.d0, q, nxi, cd, 1, 0.d0, muwk, 1)
lkhd = lkhd - ddot (nxi, cdnew, 1, muwk, 1) / 2.d0
mrs(1) = lkhd
mrs(2) = trc

return
end
