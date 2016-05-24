kron = function(mag, n, e = 0, rk = 2.5, r=1e10, re = 1){
    # kron
    bn = qgamma(0.5,2*n)
    x = bn*((r/re)^(1/n))
    krad = suppressWarnings((re/(bn^n))*(igamma(x,3*n)/igamma(x,2*n)))
    kron = sersic(mag=mag, re=re, n=n, e=e, r=rk*krad)
    return(kron)
}

