
#:::::::::::::
#   hzdaux1
#:::::::::::::

subroutine  hzdaux1 (cd, nxis, q, nxi, qdrs, nqd, qdwt, nx, mchpr, wt, v, vwk, jpvt)

integer  nxis, nxi, nqd, nx, jpvt(*)
double precision  cd(*), q(nxi,*), qdrs(nqd,nxis,*), qdwt(nqd,*), mchpr, wt(nqd,*),
                  v(nxis,*), vwk(nxis,*)

integer  i, j, k, kk, rkv
double precision  ddot

#   Initialization
for (kk=1;kk<=nx;kk=kk+1) {
    for (i=1;i<=nqd;i=i+1)
        wt(i,kk) = qdwt(i,kk) * dexp (ddot (nxis, qdrs(i,1,kk), nqd, cd, 1))
}
#   H matrix
call  dset (nxis*nxis, 0.d0, v, 1)
for (kk=1;kk<=nx;kk=kk+1) {
    for (i=1;i<=nxis;i=i+1) {
        for (j=i;j<=nxis;j=j+1) {
            vwk(i,j) = 0.d0
            for (k=1;k<=nqd;k=k+1)
                vwk(i,j) = vwk(i,j) + wt(k,kk) * qdrs(k,i,kk) * qdrs(k,j,kk)
        }
    }
    call  daxpy (nxis*nxis, 1.d0, vwk, 1, v, 1)
}
for (i=1;i<=nxi;i=i+1) {
    for (j=i;j<=nxi;j=j+1)  v(i,j) = v(i,j) + q(i,j)
}
#   Cholesky factorization
for (i=1;i<=nxis;i=i+1)  jpvt(i) = 0
call  dchdc (v, nxis, nxis, vwk, jpvt, 1, rkv)
while (v(rkv,rkv)<v(1,1)*dsqrt(mchpr))  rkv = rkv - 1
for (i=rkv+1;i<=nxis;i=i+1) {
    v(i,i) = v(1,1)
    call  dset (i-rkv-1, 0.d0, v(rkv+1,i), 1)
}

return
end


#:::::::::::::
#   hzdaux2
#:::::::::::::

subroutine  hzdaux2 (v, nxis, jpvt, r, nr, se)

double precision  v(nxis,*), r(nxis,*), se(*)
integer  nxis, jpvt(*), nr

double precision  ddot
integer  i, infowk

for (i=1;i<=nr;i=i+1) {
    call  dprmut (r(1,i), nxis, jpvt, 0)
    call  dtrsl (v, nxis, nxis, r(1,i), 11, infowk)
    se(i) = dsqrt (ddot (nxis, r(1,i), 1, r(1,i), 1))
}

return
end


#::::::::::
#   hrkl
#::::::::::

subroutine  hrkl (cd, nxis, qdrs, nqd, nx, qdwt, wt0, mchpr, wt, mu, mu0, v,
                  jpvt, wk, cdnew, wtnew, prec, maxiter, info)

integer  nxis, nqd, nx, jpvt(*), maxiter, info
double precision  cd(*), qdrs(nqd,nxis,*), qdwt(nqd,*), wt0(nqd,*), mchpr, wt(nqd,*),
                  mu(*), mu0(*), v(nxis,*), wk(*), cdnew(*), wtnew(nqd,*), prec

integer  i, j, k, kk, idamax, iter, flag, infowk
double precision  tmp, ddot, dasum, rkl, mumax, rklnew, disc, disc0


#   Initialization
info = 0
call  dset (nxis, 0.d0, mu0, 1)
for (kk=1;kk<=nx;kk=kk+1) {
    for (i=1;i<=nxis;i=i+1) 
        mu0(i) = mu0(i) + ddot (nqd, wt0(1,kk), 1, qdrs(1,i,kk), 1)
}
rkl = 0.d0
for (kk=1;kk<=nx;kk=kk+1) {
    for (i=1;i<=nqd;i=i+1) {
        tmp = ddot (nxis, qdrs(i,1,kk), nqd, cd, 1)
        wt(i,kk) = qdwt(i,kk) * dexp (tmp)
        rkl = rkl + (wt(i,kk) - wt0(i,kk)*tmp)
    }
}
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
            mu(i) = mu(i) - ddot (nqd, wt(1,kk), 1, qdrs(1,i,kk), 1)
            for (j=i;j<=nxis;j=j+1) {
                for (k=1;k<=nqd;k=k+1)
                    v(i,j) = v(i,j) + wt(k,kk) * qdrs(k,i,kk) * qdrs(k,j,kk)
            }
        }
    }
    call  daxpy (nxis, 1.d0, mu0, 1, mu, 1)
    mumax = dabs(mu(idamax(nxis, mu, 1)))
    #   Cholesky factorization
    for (i=1;i<=nxis;i=i+1)  jpvt(i) = 0
    call  dmcdc (v, nxis, nxis, wk, jpvt, infowk)
    #   Update coefficients
    repeat {
        call  dcopy (nxis, mu, 1, cdnew, 1)
        call  dprmut (cdnew, nxis, jpvt, 0)
        call  dtrsl (v, nxis, nxis, cdnew, 11, infowk)
        call  dtrsl (v, nxis, nxis, cdnew, 01, infowk)
        call  dprmut (cdnew, nxis, jpvt, 1)
        call  daxpy (nxis, 1.d0, cd, 1, cdnew, 1)
        rklnew = 0.d0
        for (kk=1;kk<=nx;kk=kk+1) {
            for (i=1;i<=nqd;i=i+1) {
                tmp = ddot (nxis, qdrs(i,1,kk), nqd, cdnew, 1)
                if (tmp>3.d2) {
                    flag = flag + 1
                    break
                }
                wtnew(i,kk) = qdwt(i,kk) * dexp (tmp)
                rklnew = rklnew + (wtnew(i,kk) - wt0(i,kk)*tmp)
            }
        }
        if (flag==1) {
            #   Reset iteration with uniform starting value
            call  dset (nxis, 0.d0, cd, 1)
            rkl = dasum (nqd*nx, qdwt, 1)
            call  dcopy (nqd*nx, qdwt, 1, wt, 1)
            iter = 0
            break
        }
        if (rklnew-rkl<1.d1*(1.d0+dabs(rkl))*mchpr)  break
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
    disc = dmax1 (disc, (mumax/(1.d0+rkl))**2)
    disc0 = dmax1 ((mumax/(1.d0+rkl))**2, dabs(rkl-rklnew)/(1+dabs(rkl)))
    #   Set to new values
    call  dcopy (nxis, cdnew, 1, cd, 1)
    call  dcopy (nqd*nx, wtnew, 1, wt, 1)
    rkl = rklnew
    #   Check convergence
    if (disc0<prec)  break
    if (disc<prec)  break
    if (iter<maxiter)  next
    if (flag==0) {
        #   Reset iteration with uniform starting value
        call  dset (nxis, 0.d0, cd, 1)
        rkl = dasum (nqd*nx, qdwt, 1)
        call  dcopy (nqd*nx, qdwt, 1, wt, 1)
        iter = 0
        flag = 2
    }
    else {
        info = 2
        break
    }
}
#   return projection density on the mesh
for (kk=1;kk<=nx;kk=kk+1) {
    for (i=1;i<=nqd;i=i+1)
        wt0(i,kk) = dexp (ddot (nxis, qdrs(i,1,kk), nqd, cd, 1))
}

return
end


#:::::::::::::
#   coxaux
#:::::::::::::

subroutine  coxaux (cd, nn, q, nxiz, qdrs, nqd, nt, twt, mchpr, qdwt,
                    wt, wtsum, muwk, v, vwk, jpvt)

integer  nn, nxiz, nqd, nt, jpvt(*)
double precision  cd(*), q(nxiz,*), qdrs(nqd,*), twt(*), mchpr, qdwt(nqd,*),
                  wt(nqd,*), wtsum(*), muwk(*), v(nn,*), vwk(nn,*)

integer  i, j, k, m, rkv
double precision  ddot, tmp

#   Initialization
call  dset(nt, 0.d0, wtsum, 1)
for (i=1;i<=nqd;i=i+1) {
    tmp = dexp (ddot (nn, qdrs(i,1), nqd, cd, 1))
    for (m=1;m<=nt;m=m+1) {
        wt(i,m) = qdwt(i,m) * tmp
        wtsum(m) = wtsum(m) + wt(i,m)
    }
}
#   H matrix
call  dset(nn*nn, 0.d0, v, 1)
for (m=1;m<=nt;m=m+1) {
    for (i=1;i<=nn;i=i+1)
        muwk(i) = ddot (nqd, wt(1,m), 1, qdrs(1,i), 1) / wtsum(m)
    for (i=1;i<=nn;i=i+1) {
        for (j=i;j<=nn;j=j+1) {
            vwk(i,j) = 0.d0
            for (k=1;k<=nqd;k=k+1)
                vwk(i,j) = vwk(i,j) + wt(k,m) * qdrs(k,i) * qdrs(k,j)
            vwk(i,j) = vwk(i,j) / wtsum(m) - muwk(i) * muwk(j)
        }
    }
    call  daxpy (nn*nn, twt(m), vwk, 1, v, 1)
}
for (i=1;i<=nxiz;i=i+1) {
    for (j=i;j<=nxiz;j=j+1)  v(i,j) = v(i,j) + q(i,j)
}
#   Cholesky factorization
for (i=1;i<=nn;i=i+1)  jpvt(i) = 0
call  dchdc (v, nn, nn, vwk, jpvt, 1, rkv)
while (v(rkv,rkv)<v(1,1)*dsqrt(mchpr))  rkv = rkv - 1
for (i=rkv+1;i<=nn;i=i+1) {
    v(i,i) = v(1,1)
    call  dset (i-rkv-1, 0.d0, v(rkv+1,i), 1)
}

return
end
