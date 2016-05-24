petroindex = function(r, n, re = 1){
    # petrosian index
    bn = qgamma(0.5,2*n)
    x = bn*((r/re)^(1/n))
    i = (2*n*igamma(x,2*n))/((exp(1)^(-x))*(x^(2*n)))
    return(i)
}

