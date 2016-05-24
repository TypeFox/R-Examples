
#::::::::::
#   drkl
#::::::::::

subroutine  drkl (cd, nxis, qdrs, nqd, nt, bwt, qdwt, wt0, mchpr, prec, maxiter,
                  jpvt, wk, info)

integer  nxis, nqd, nt, maxiter, jpvt(*), info
double precision  cd(*), qdrs(nqd,*), bwt(*), qdwt(nt,*), wt0(*), mchpr, prec, wk(*)

integer  imrs, iwt, iwtsum, imu, imuwk, iv, ivwk, icdnew, iwtnew, iwtnewsum

imrs = 1
iwt = imrs + nxis
iwtsum = iwt + nt*nqd
imu = iwtsum + nt
imuwk = imu + nxis
iv = imuwk + nxis
ivwk = iv + nxis*nxis
icdnew = ivwk + nxis*nxis
iwtnew = icdnew + nxis
iwtnewsum = iwtnew + nt*nqd

call  drkl1 (cd, nxis, qdrs, nqd, nt, bwt, qdwt, wt0, mchpr, prec, maxiter,
             wk(imrs), wk(iwt), wk(iwtsum), wk(imu), wk(imuwk), wk(iv), wk(ivwk),
             jpvt, wk(icdnew), wk(iwtnew), wk(iwtnewsum), info)

return
end


#::::::::::
#   drkl1
#::::::::::

subroutine  drkl1 (cd, nxis, qdrs, nqd, nt, bwt, qdwt, wt0, mchpr, prec, maxiter,
                   mrs, wt, wtsum, mu, muwk, v, vwk, jpvt, cdnew, wtnew, wtnewsum,
                   info)

integer  nxis, nqd, nt, maxiter, jpvt(*), info
double precision  cd(*), qdrs(nqd,*), bwt(*), qdwt(nt,*), wt0(*), mchpr, mrs(*),
                  wt(nt,*), wtsum(*), mu(*), muwk(*), v(nxis,*), vwk(nxis,*),
                  cdnew(*), wtnew(nt,*), wtnewsum(*), prec

integer  i, j, k, m, iter, flag, idamax, infowk
double precision  tmp, ddot, rkl, rklnew, mumax, disc, disc0

#   Initialization
info = 0
call  dset (nxis, 0.d0, mrs, 1)
for (m=1;m<=nt;m=m+1) {
    for (i=1;i<=nqd;i=i+1)  wt(m,i) = qdwt(m,i) * wt0(i)
    for (i=1;i<=nxis;i=i+1)  muwk(i) = ddot (nqd, qdrs(1,i), 1, wt(m,1), nt)
    call  daxpy (nxis, bwt(m), muwk, 1, mrs, 1)
    wtsum(m) = 0.d0
}
for (i=1;i<=nqd;i=i+1) {
    tmp = dexp (ddot (nxis, qdrs(i,1), nqd, cd, 1))
    for (m=1;m<=nt;m=m+1) {
        wt(m,i) = qdwt(m,i) * tmp
        wtsum(m) = wtsum(m) + wt(m,i)
    }
}
rkl = 0.d0
for (m=1;m<=nt;m=m+1) {
    tmp = 0.d0
    for (i=1;i<=nqd;i=i+1) {
        disc = wt0(i) * qdwt(m,i)
        tmp = tmp + dlog (disc/wt(m,i)) * disc
    }
    rkl = rkl + bwt(m) * (tmp + dlog (wtsum(m)))
}
iter = 0
flag = 0
#   Iteration
repeat {
    iter = iter + 1
    # Calculate hessian and gradient
    call  dset(nxis, 0.d0, mu, 1)
    call  dset(nxis*nxis, 0.d0, v, 1)
    for (m=1;m<=nt;m=m+1) {
        for (i=1;i<=nxis;i=i+1)
            muwk(i) = - ddot (nqd, wt(m,1), nt, qdrs(1,i), 1) / wtsum(m)
        for (i=1;i<=nxis;i=i+1) {
            for (j=i;j<=nxis;j=j+1) {
                vwk(i,j) = 0.d0
                for (k=1;k<=nqd;k=k+1)
                    vwk(i,j) = vwk(i,j) + wt(m,k) * qdrs(k,i) * qdrs(k,j)
                vwk(i,j) = vwk(i,j) / wtsum(m) - muwk(i) * muwk(j)
            }
        }
        call  daxpy (nxis, bwt(m), muwk, 1, mu, 1)
        call  daxpy (nxis*nxis, bwt(m), vwk, 1, v, 1)
    }
    call  daxpy (nxis, 1.d0, mrs, 1, mu, 1)
    mumax = dabs(mu(idamax(nxis, mu, 1)))
    #   Cholesky factorization
    for (i=1;i<=nxis;i=i+1)  jpvt(i) = 0
    call  dmcdc (v, nxis, nxis, muwk, jpvt, infowk)
    #   Update coefficients
    repeat {
        call  dcopy (nxis, mu, 1, cdnew, 1)
        call  dprmut (cdnew, nxis, jpvt, 0)
        call  dtrsl (v, nxis, nxis, cdnew, 11, infowk)
        call  dtrsl (v, nxis, nxis, cdnew, 01, infowk)
        call  dprmut (cdnew, nxis, jpvt, 1)
        call  daxpy (nxis, 1.d0, cd, 1, cdnew, 1)
        call  dset (nt, 0.d0, wtnewsum, 1)
        for (i=1;i<=nqd;i=i+1) {
            tmp = ddot (nxis, qdrs(i,1), nqd, cdnew, 1)
            if (tmp>3.d2) {
                flag = flag + 1
                break
            }
            for (m=1;m<=nt;m=m+1) {
                wtnew(m,i) = qdwt(m,i) * dexp (tmp)
                wtnewsum(m) = wtnewsum(m) + wtnew(m,i)
            }
        }
        rklnew = 0.d0
        for (m=1;m<=nt;m=m+1) {
            tmp = 0.d0
            for (i=1;i<=nqd;i=i+1) {
                disc = wt0(i) * qdwt(m,i)
                tmp = tmp + dlog (disc/wtnew(m,i)) * disc
            }
            rklnew = rklnew + bwt(m) * (tmp + dlog (wtnewsum(m)))
        }
        if (flag==1) {
            #   Reset iteration with uniform starting value
            call  dset (nxis, 0.d0, cd, 1)
            call  dcopy (nt*nqd, qdwt, 1, wt, 1)
            call  dset (nt, 0.d0, wtsum, 1)
            for (m=1;m<=nt;m=m+1) {
                for (i=1;i<=nqd;i=i+1)  wtsum(m) = wtsum(m) + wt(m,i)
            }
            rkl = 0.d0
            for (m=1;m<=nt;m=m+1) {
                tmp = 0.d0
                for (i=1;i<=nqd;i=i+1)
                    tmp = tmp + dlog (wt0(i)) * wt0(i) *  qdwt(m,i)
                rkl = rkl + bwt(m) * (tmp + dlog (wtsum(m)))
            }
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
    for (i=1;i<=nqd;i=i+1){
        for (m=1;m<=nt;m=m+1)
            disc = dmax1 (disc, dabs(wt(m,i)-wtnew(m,i))/(1.d0+dabs(wt(m,i))))
    }
    disc = dmax1 (disc, (mumax/(1.d0+rkl))**2)
    disc = dmax1 (disc, dabs(rkl-rklnew)/(1.d0+dabs(rkl)))
    #   Check convergence
    if (disc<prec)  break
    #   Set to new values
    call  dcopy (nxis, cdnew, 1, cd, 1)
    call  dcopy (nt*nqd, wtnew, 1, wt, 1)
    call  dcopy (nt, wtnewsum, 1, wtsum, 1)
    rkl = rklnew
    if (iter<maxiter)  next
    if (flag==0) {
        #   Reset iteration with uniform starting value
        call  dset (nxis, 0.d0, cd, 1)
        call  dcopy (nt*nqd, qdwt, 1, wt, 1)
        call  dset (nt, 0.d0, wtsum, 1)
        for (m=1;m<=nt;m=m+1) {
            for (i=1;i<=nqd;i=i+1)  wtsum(m) = wtsum(m) + wt(m,i)
        }
        rkl = 0.d0
        for (m=1;m<=nt;m=m+1) {
            tmp = 0.d0
            for (i=1;i<=nqd;i=i+1)
                tmp = tmp + dlog (wt0(i)) * wt0(i) *  qdwt(m,i)
            rkl = rkl + bwt(m) * (tmp + dlog (wtsum(m)))
        }
        iter = 0
        flag = 2
    }
    else {
        info = 2
        break
    }
}
#   return projection density on the mesh
for (i=1;i<=nqd;i=i+1)  wt0(i) = dexp (ddot (nxis, qdrs(i,1), nqd, cd, 1))

return
end
