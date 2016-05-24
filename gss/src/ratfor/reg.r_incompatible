
#:::::::::
#   reg
#:::::::::

subroutine  reg (sr, nobs, nnull, q, nxi, y, method, alpha, varht,
                 score, dc, mchpr, v, mu, jpvt, wk, rkv, info)

double precision  sr(nobs,*), q(nxi,*), y(*), alpha, varht, score, dc(*),
                  mchpr, v(nnull+nxi,*), mu(*), wk(*)
integer  nobs, nnull, nxi, method, jpvt(*), rkv, info

double precision  ddot, dasum, rss, trc, dum
integer  i, j, nn, idamax, infowk

info = 0
nn = nnull + nxi

#   form linear system
for (i=1;i<=nn;i=i+1) {
    mu(i) = ddot (nobs, sr(1,i), 1, y, 1)
    for (j=i;j<=nn;j=j+1) {
        v(i,j) = ddot (nobs, sr(1,i), 1, sr(1,j), 1)
        if (i>nnull)  v(i,j) = v(i,j) + q(i-nnull,j-nnull)
    }
}
#   Cholesky decomposition
infowk = 0
for (i=1;i<=nn;i=i+1)  infowk = infowk + jpvt(i)
call  dchdc (v, nn, nn, wk, jpvt, 1, rkv)
j = idamax (rkv-infowk, v(infowk+1,infowk+1), nn+1)
while (v(rkv,rkv)<v(infowk+j,infowk+j)*dsqrt(mchpr))  rkv = rkv - 1
for (i=rkv+1;i<=nn;i=i+1) {
    v(i,i) = v(j,j)
    call  dset (i-rkv-1, 0.d0, v(rkv+1,i), 1)
}
#   obtain coefficients
call  dcopy (nn, mu, 1, dc, 1)
call  dprmut (dc, nn, jpvt, 0)
call  dtrsl (v, nn, nn, dc, 11, infowk)
call  dset (nn-rkv, 0.d0, dc(rkv+1), 1)
call  dtrsl (v, nn, nn, dc, 01, infowk)
call  dprmut (dc, nn, jpvt, 1)
#   early return
if (method==4)  return
#   calculate residuals
for (i=1;i<=nobs;i=i+1)  wk(i) = y(i) - ddot (nn, sr(i,1), nobs, dc, 1)
#   diagonal of smoothing matrix
if (method==5) {
    wk(nobs+1) = ddot (nobs, wk, 1, wk, 1) / dfloat (nobs)
    for (i=1;i<=nobs;i=i+1) {
        call  dcopy (nn, sr(i,1), nobs, mu, 1)
        call  dprmut (mu, nn, jpvt, 0)
        call  dtrsl (v, nn, nn, mu, 11, infowk)
        wk(i) = ddot (nn, mu, 1, mu, 1)
    }
    return    
}
if (method==3) {
    #   GML
    rss = ddot (nobs, y, 1, wk, 1)
    #   determinant
    if (nnull>0) {
        call  dqrdc (sr, nobs, nobs, nnull, wk, dum, dum, 0)
        for (i=1;i<=nxi;i=i+1) {
            call  dqrsl (sr, nobs, nobs, nnull, wk, sr(1,nnull+i),
                         dum, sr(1,nnull+i), dum, dum, dum, 01000, infowk)
        }
    }
    call  dcopy (nxi, q, nxi+1, wk, 1)
    for (i=1;i<=nxi;i=i+1) {
        for (j=i;j<=nxi;j=j+1)
            q(i,j) = q(i,j) + ddot (nobs-nnull, sr(nnull+1,nnull+i), 1,
                                    sr(nnull+1,nnull+j), 1)
    }
    for (i=1;i<=nxi;i=i+1) {
        for (j=i;j<=nxi;j=j+1) {
            sr(i,j) = q(i,j)
            sr(j,i) = q(i,j)
            q(i,j) = q(j,i)
        }
    }
    call  dcopy (nxi, wk, 1, q, nxi+1)
#    call  rs (nobs, nxi, sr, mu, 0, dum, wk, y, info)
    call  dsyev ('n', 'u', nxi, sr, nobs, mu, wk, 3*nxi, info)
    trc = 0.d0
    for (i=1;i<=rkv-nnull;i=i+1)  trc = trc + dlog (mu(nxi-i+1))
#    call  rs (nxi, nxi, q, mu, 0, dum, wk, y, info)
    call  dsyev ('n', 'u', nxi, q, nxi, mu, wk, 3*nxi, info)
    for (i=1;i<=rkv-nnull;i=i+1)  trc = trc - dlog (mu(nxi-i+1))
    #   return values
    score = rss / dfloat (nobs) * dexp (trc/dfloat(nobs-nnull))
    varht = rss / dfloat (nobs-nnull)
}
else {
    #   GCV or Cp
    rss = ddot (nobs, wk, 1, wk, 1) / dfloat (nobs)
    #   trace
    for (i=1;i<=nobs;i=i+1) {
        call  dcopy (nn, sr(i,1), nobs, mu, 1)
        call  dprmut (mu, nn, jpvt, 0)
        call  dtrsl (v, nn, nn, mu, 11, infowk)
        wk(i) = ddot (nn, mu, 1, mu, 1)
    }
    trc = dasum (nobs, wk, 1) / dfloat (nobs)
    #   return values
    if (method==2) {
        score = rss / (1.d0-alpha*trc)**2
        varht = rss / (1.d0-trc)
    }
    else  score = rss + 2.d0 * varht * alpha * trc
}
wk(1) = rss
wk(2) = trc

return
end


#::::::::::::
#   regaux
#::::::::::::

subroutine  regaux (v, nn, jpvt, rkv, r, nr, sms, nnull, wk)

double precision  v(nn,*), r(nn,*), sms(nnull,*), wk(nn,*)
integer  nn, jpvt(*), rkv, nr, nnull

double precision  ddot
integer  i, j, infowk

#   drcr
for (i=1;i<=nr;i=i+1) {
    call  dprmut (r(1,i), nn, jpvt, 0)
    call  dtrsl (v, nn, nn, r(1,i), 11, infowk)
    if (nn-rkv>0)  call  dset (nn-rkv, 0.d0, r(rkv+1,i), 1)
    call  dtrsl (v, nn, nn, r(1,i), 01, infowk)
    call  dprmut (r(1,i), nn, jpvt, 1)
}
#   sms
call  dset (nn*nnull, 0.d0, wk, 1)
call  dset (nnull, 1.d0, wk, nn+1)
for (i=1;i<=nnull;i=i+1)
    call  dtrsl (v, nn, nn, wk(1,i), 11, infowk)
for (i=1;i<=nnull;i=i+1) {
    for (j=i;j<=nnull;j=j+1) {
        sms(i,j) = ddot (nn, wk(1,i), 1, wk(1,j), 1)
        sms(j,i) = sms(i,j)
    }
}

return
end
