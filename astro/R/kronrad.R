kronrad = function(n, r = 1e10, re = 1){
    # kron radii
    bn = qgamma(0.5,2*n)
    x = bn*((r/re)^(1/n))
    krad = suppressWarnings((re/(bn^n))*(igamma(x,3*n)/igamma(x,2*n)))
    return(krad)
}

